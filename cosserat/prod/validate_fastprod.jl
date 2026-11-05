## Check fastprod.jl against run_simulation_plastic_remeshed (../src/rod_contact2.jl)
#
#   julia validate_fastprod.jl        (or sbatch validate_fastprod.sh)
#
# Full pipeline on V6 data (walls, self-contact, remeshing, data-driven clamped ends,
# growth, remodeling) for a short simulation time T: the growth rate is rescaled by 1/T,
# so the whole measured elongation, and the remeshing it triggers, still happens.
# Compares the saved xsaves / Qsaves of both kernels and times them.

import Pkg; Pkg.activate("../src/rods")
include("../src/rod_contact2.jl")
using LinearAlgebra
using Random
using MAT
using Interpolations
using Printf
include("../src/fastkernel.jl")
include("fastprod.jl")
include("prod_common.jl")

function prod_bcs(smp, T, L) # clamped, data-driven ends for run_simulation_plastic_remeshed
    x_itp = smp.x_itp
    Q1_itp = smp.Q1_itp
    Qend_itp = smp.Qend_itp
    function bcs(t, x, Q)
        xp = copy(x)
        Qp = copy(Q)
        Qp[1, :, :] = copy(Q1_itp(t / T))
        Qp[end, :, :] = copy(Qend_itp(t / T))
        xp[1, :] = x_itp(t / T)[1, :] / L
        xp[end, :] = x_itp(t / T)[end, :] / L
        return xp, Qp
    end
    return bcs
end

function ext_forces(x, neighbour_walls) # wall forces for run_simulation_plastic_remeshed
    f_ext = zeros(size(x, 1), 3)
    N = size(x, 1) - 1
    for i = 2:N
        force_orient, dist = get_force_components(x[i, :], neighbour_walls[i])
        f_ext[i, :] = 10000 * dist * (1 + 00 * dist) * force_orient
    end
    return f_ext
end
ext_torques(t::Float64, x::Matrix{Float64}, Q::Array{Float64,3}) = zeros(size(Q, 1), 3)

ok = true
for (sample_id, T, tw, p) in [(5, 0.1, 1.0, 10.0), (2, 0.1, 0.3, 0.0)]
    smp = setup_sample(sample_id)
    r = prepare_run(smp, T, tw, p)
    n_saves = 20

    rod = make_Rod(r.x_start, r.S, r.St, r.B, r.Bt, r.gamma_p, r.gamma_t)
    t0 = time()
    xp, Qp = run_simulation_plastic_remeshed(copy(r.x_start), copy(r.Qstart), r.dt, r.Ns, T, n_saves,
        rod, r.S, r.St, prod_bcs(smp, T, r.L), r.circum_radii, r.alpha, r.walls, r.cutoff,
        ext_forces, ext_torques, r.growth_rate, p, r.dxmax)
    tp = time() - t0

    t0 = time()
    xf, Qf = run_fastprod(copy(r.x_start), copy(r.Qstart), r.dt, r.Ns, T, n_saves,
        r.S, r.St, r.B, r.Bt, r.gamma_p, r.gamma_t, make_bc_fast(smp, T, r.L), r.circum_radii,
        r.alpha, r.walls, r.cutoff, r.W1, r.growth_rate, p, r.dxmax; material_u0=true)
    tf = time() - t0

    steps = (n_saves - 1) * round(Int, T / r.dt / n_saves) + 1
    @printf("sample %d (%s), T=%.2f tw=%.1f p=%.0f, %d steps:\n", sample_id, smp.cfg.ref, T, tw, p, steps)
    for k in [2, n_saves ÷ 2, n_saves]
        ex = maximum(abs.(xf[k, :, :] .- xp[k, :, :]))
        eQ = maximum(abs.(Qf[k, :, :, :] .- Qp[k, :, :, :]))
        @printf("    save %2d: max |Δx| %.2e (rod length 1), max |ΔQ| %.2e\n", k, ex, eQ)
    end
    ex = maximum(abs.(xf .- xp)); eQ = maximum(abs.(Qf .- Qp))
    @printf("    all saves: max |Δx| %.2e, max |ΔQ| %.2e\n", ex, eQ)
    @printf("    time: prod %.1f s (%.0f μs/step), fast %.1f s (%.0f μs/step), %.0f×\n",
            tp, 1e6 * tp / steps, tf, 1e6 * tf / steps, tp / tf)
    global ok &= (ex < 1e-8 && eQ < 1e-8)
end
println(ok ? "~fastprod matches prod~" : "!! MISMATCH between fastprod and prod !!")
