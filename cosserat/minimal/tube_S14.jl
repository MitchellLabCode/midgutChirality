## Growth of the S14 rods in a hard-walled tube or without confinement (minimal version of
## noah/rods/reduced/tube_full.jl, restricted to the S14 family, hard / free, uniform growth)
#
# Overdamped Cosserat rod (RK4), clamped at x = 0 and x = Le ez with fixed end frames, grown
# uniformly from L0 = 1.05 Le to L_end, with Maxwell remodeling of the rest curvature at the
# fluidity η, self-contact (circumradius r_contact), remeshing, and rotational noise.
#   hard : centerline confined to ρ ≤ R by a stiff wall force -W (ρ - R) ρ̂ per unit length
#   free : no confinement
# Initial shapes: the lowest clamped pair, ψ = A (e1 + 0.5 i e2), and its mirror image,
# interpolated by the initial asymmetry α ∈ [-1, 1] (α = 0 planar); zero link.
# Lengths in units of Le, times in τE (code units B1 = Le = γ = 1).
#
# Usage (from this folder):
#   julia -p 7 tube_S14.jl hard            # or free, or both (default)
# Options (environment variables):
#   TUBE_LEND   final length L_end/Le                     (default 1.5)
#   TUBE_NALPHA number of α values in [-1, 1]              (default 7)
#   TUBE_NSAVE  saved shapes resampled on this many edges  (default 120)
#   TUBE_SIGMA  rotational noise strength (0: no noise)    (default 1e-3)
#   TUBE_ETAS   fluidities η τE                           (default "0 10 100 1000")
#   TUBE_ONLY_ETA  run only this η τE (one of TUBE_ETAS; the seeds depend on its position)
#   TUBE_R      tube radius R/Le                           (default 0.3)
#   TUBE_DATE   date stamp of the output folder            (default 20261020)
# Output: ./<date>_tube_S14_<hard|free>_L<L_end>[_na<n>][_nonoise]/eta<η>_a<i>.mat (resumable),
# with xsaves (200 × (Ns+1) × 3), Qsaves, t, Lrest, N, Wr_fuller, Tw, rho_max and the parameters.

using Distributed

@everywhere begin import Pkg; Pkg.activate("../src/rods")
    Pkg.instantiate();
end
@everywhere begin
    include("../src/rod_contact2.jl")   # interpolate_rodstate, Plane (used by fastprod.jl)
    using LinearAlgebra
    using Random
    using MAT
    include("../src/fastkernel.jl")     # FastRod, rk4_step!, static-array helpers
    include("../prod/fastprod.jl")      # remesh_fast, neighbour_points!, add_contact_forces!
end

@everywhere begin
    # ---- parameters ------------------------------------------------------------
    const Le            = 1.0
    const L_start       = 1.05
    const L_end_tube    = parse(Float64, get(ENV, "TUBE_LEND", "1.5"))
    const g_rate        = 1.0          # g τE
    const B1            = 1.0
    const B3            = 1.0          # B3/B1 = 1
    const S_stiff       = 6e3          # stretch/shear stiffness: axial strain ≪ excess length
    const gamma_R       = 1.0          # rotational drag per unit length
    const gamma_L       = 1.0          # translational drag per unit length
    const sigma         = parse(Float64, get(ENV, "TUBE_SIGMA", "1e-3"))   # noise strength
    const material_u0   = true         # rest curvature stored in the material frame
    const N0            = 40           # initial number of edges
    const n_saves       = 200
    const dt_safety     = 1.0          # dt = dt_safety / (fastest relaxation rate)
    const tube_date     = get(ENV, "TUBE_DATE", "20261020")
    const tube_radius   = parse(Float64, get(ENV, "TUBE_R", "0.3"))   # hard-wall radius for the centerline
    const W_wall        = 1e4          # wall stiffness per unit length
    const r_contact     = 0.14         # self-contact (circum)radius
    const alpha_contact = 1e4          # contact stiffness per unit length
    # η τE values; the seed of a run depends on the position of its η in this list
    const tube_etas     = parse.(Float64, split(get(ENV, "TUBE_ETAS", "0 10 100 1000")))
    const tube_nalpha   = parse(Int, get(ENV, "TUBE_NALPHA", "7"))
    const tube_alphas   = collect(range(-1.0, 1.0, length=tube_nalpha))
    const tube_confs    = ["hard", "soft", "free"]   # (soft: not in this version; kept for the seeds)
    const Ns_save       = parse(Int, get(ENV, "TUBE_NSAVE", "120"))
    const ez            = [0.0, 0.0, 1.0]

    # ---- geometry helpers --------------------------------------------------------
    function rotmat(theta::Float64, u::Vector{Float64})
        """rotation by theta about the unit vector u"""
        U = [0.0 -u[3] u[2]; u[3] 0.0 -u[1]; -u[2] u[1] 0.0]
        return I + sin(theta) * U + (1 - cos(theta)) * U * U
    end

    function minrot(a::Vector{Float64}, b::Vector{Float64})
        """minimal rotation sending the unit vector a onto the unit vector b"""
        c = cross(a, b)
        s = norm(c)
        s < 1e-14 && return Matrix{Float64}(I, 3, 3)
        return rotmat(atan(s, dot(a, b)), c / s)
    end

    function clamped_modes(s::Vector{Float64}, w::Vector{Float64}, L::Float64, nmodes::Int)
        """lowest `nmodes` eigenfunctions of -∂s² with e(0) = e(L) = 0 and ∫e = 0 at the points s
        (quadrature weights w), normalized by (2/L)∫e² = 1, sorted by wavenumber: sin(2πn s/L) and
        (cos k x - cos ϑ)/sin ϑ with tan ϑ = ϑ, ϑ = kL/2, x = s - L/2"""
        ks = Float64[]; kinds = Symbol[]; thetas = Float64[]
        for n = 1:nmodes
            push!(ks, 2π * n / L); push!(kinds, :sin); push!(thetas, 0.0)
        end
        for m = 1:nmodes
            θ = (m + 0.5) * π - 1 / ((m + 0.5) * π)
            for _ = 1:50
                θ -= (tan(θ) - θ) / tan(θ)^2
            end
            push!(ks, 2θ / L); push!(kinds, :cos); push!(thetas, θ)
        end
        order = sortperm(ks)[1:nmodes]
        E = zeros(length(s), nmodes)
        x = s .- L / 2
        for (j, o) in enumerate(order)
            if kinds[o] == :sin
                e = sin.(ks[o] .* s)
            else
                e = (cos.(ks[o] .* x) .- cos(thetas[o])) ./ sin(thetas[o])
            end
            e ./= sqrt(2 / L * sum(w .* e .^ 2))
            E[:, j] = e
        end
        return E
    end

    # ---- diagnostics ---------------------------------------------------------------
    function edge_tangents(x::Matrix{Float64})
        d = x[2:end, :] - x[1:end-1, :]
        ls = sqrt.(sum(abs2, d, dims=2))
        return d ./ ls, vec(ls)
    end

    function writhe_fuller(x::Matrix{Float64})
        """Fuller's formula relative to the axis ez: Wr = (1/2π) ∫ ds (ez × t)·t' / (1 + ez·t)"""
        t, _ = edge_tangents(x)
        wr = 0.0
        for k = 1:size(t, 1)-1
            tb = normalize(t[k, :] + t[k+1, :])
            wr += dot(cross(ez, tb), t[k+1, :] - t[k, :]) / (1 + tb[3])
        end
        return wr / (2π)
    end

    function twist_total(Q::Array{Float64,3})
        """Tw = (1/2π) Σ twist angle between successive frames, after parallel transport"""
        tw = 0.0
        for k = 1:size(Q, 1)-1
            t1 = Q[k, 3, :]; t2 = Q[k+1, 3, :]
            d1t = minrot(t1, t2) * Q[k, 1, :]
            d1n = Q[k+1, 1, :]
            tw += atan(dot(cross(d1t, d1n), t2), dot(d1t, d1n))
        end
        return tw / (2π)
    end

    # ---- initial shape -------------------------------------------------------------
    function initial_state(alpha::Float64)
        """S14 shape for asymmetry α: ψ = tx + i ty = A Σ a_j e_j with a = (1, 0.5 i) at α = -1 and
        its mirror image a* at α = +1, a_α = (1-α)/2 a + (1+α)/2 a*; A from the exact end shortening
        Σ ℓ sqrt(1 - A²|ψ|²) = Le; frames by parallel transport plus a uniform twist (zero link)"""
        N = N0
        ell = L_start / N
        s = collect(((1:N) .- 0.5) .* ell)          # edge midpoints
        w = fill(ell, N)
        E = clamped_modes(s, w, L_start, 2)
        a = ComplexF64[1.0, im * 0.5]
        a_alpha = (1 - alpha) / 2 .* a .+ (1 + alpha) / 2 .* conj.(a)
        psi = E * a_alpha
        psi .-= sum(psi) / N                        # exact discrete ∫ψ = 0: X(L) = 0

        axial(A) = ell * sum(sqrt.(1 .- A^2 .* abs2.(psi)))
        A_lo, A_hi = 0.0, 0.999 / maximum(abs.(psi))
        axial(A_hi) < Le || error("initial shape cannot absorb the excess length")
        for _ = 1:200
            A_mid = 0.5 * (A_lo + A_hi)
            axial(A_mid) > Le ? (A_lo = A_mid) : (A_hi = A_mid)
        end
        A = 0.5 * (A_lo + A_hi)
        t = hcat(A .* real.(psi), A .* imag.(psi), sqrt.(1 .- A^2 .* abs2.(psi)))

        x = zeros(N + 1, 3)
        for k = 1:N
            x[k+1, :] = x[k, :] + ell * t[k, :]
        end
        x[end, :] = [0.0, 0.0, Le]

        F = zeros(N, 3, 3)                          # columns = directors
        F[1, :, :] = minrot(ez, t[1, :])
        for k = 1:N-1
            F[k+1, :, :] = minrot(t[k, :], t[k+1, :]) * F[k, :, :]
        end
        G = minrot(ez, t[N, :])
        tN = t[N, :]
        delta = atan(dot(cross(F[N, :, 1], G[:, 1]), tN), dot(F[N, :, 1], G[:, 1]))
        Q = zeros(N, 3, 3)                          # rows = directors (kernel convention)
        for k = 1:N
            Fk = rotmat(delta * (k - 1) / (N - 1), t[k, :]) * F[k, :, :]
            Q[k, :, :] = Fk'
        end
        return x, Q
    end

    # ---- hard wall -----------------------------------------------------------------
    function hard_wall!(F, x, R, W)
        """force per node -W (ρ - R) ρ̂ × (Voronoi length) on interior nodes beyond ρ = R"""
        @inbounds for i = 2:length(x)-1
            rho = sqrt(x[i][1]^2 + x[i][2]^2)
            if rho > R
                lv = 0.5 * (norm(x[i] - x[i-1]) + norm(x[i+1] - x[i]))
                c = -W * (rho - R) * lv / rho
                F[i] += V3(c * x[i][1], c * x[i][2], 0.0)
            end
        end
        return nothing
    end

    # ---- one run -------------------------------------------------------------------
    function run_tube(conf::String, eta::Float64, ia::Int)
        alpha = tube_alphas[ia]
        ctag = (conf == "hard" && tube_radius != 0.3) ? "hard_R$(tube_radius)" : conf
        ltag = (L_end_tube == 2.5 ? "" : "_L$(L_end_tube)") *
               (tube_nalpha == 7 ? "" : "_na$(tube_nalpha)") *
               (sigma == 1e-3 ? "" : sigma == 0 ? "_nonoise" : "_sig$(sigma)")
        folder_name = "./$(tube_date)_tube_S14_$(ctag)$(ltag)"
        out_filename = folder_name * "/eta$(eta)_a$(ia).mat"
        isfile(out_filename) && return out_filename
        mkpath(folder_name)

        seed = 1100000 + 10000 * findfirst(==(conf), tube_confs) + 100 * findfirst(==(eta), tube_etas) + ia
        rng = MersenneTwister(seed)

        x0, Q0 = initial_state(alpha)
        x = to_static_x(x0)
        Q = to_static_Q(Q0)
        N = length(Q)
        Q1 = Q[1]; QN = Q[N]
        x_end = V3(0.0, 0.0, Le)
        function bc!(x, Q)
            x[1] = zero(V3); x[end] = x_end
            Q[1] = Q1; Q[end] = QN
            return nothing
        end

        ell0 = L_start / N
        dxmax = 2.0 * ell0                      # remesh threshold
        l0s = fill(ell0, N)
        u0 = zeros(V3, N - 1)                   # straight, untwisted rest state
        w = ones(N)
        rod = FastRod(N, S_stiff, S_stiff, B1, B3, gamma_L * ell0, gamma_R * ell0; material_u0=material_u0)
        buf = FastBuffers(N)

        nl = Vector{Vector{Int}}(); arclen = Float64[]; fc = V3[]
        cradii = fill(r_contact, N + 1)
        alpha_node = Ref(alpha_contact * ell0)  # contact stiffness per node
        function extf!(F, xx)
            conf == "hard" && hard_wall!(F, xx, tube_radius, W_wall)
            add_contact_forces!(F, xx, nl, cradii, alpha_node[], fc)
        end
        neighbour_points!(nl, arclen, x, cradii)

        # edges never get shorter than ell0 (split at 2 ell0), so dt is set at the start
        rate = max(4 * S_stiff / (gamma_L * ell0^2), 16 * B1 / (gamma_L * ell0^4), W_wall / gamma_L,
                   alpha_contact / gamma_L)
        dt0 = dt_safety / rate
        T = log(L_end_tube / L_start) / g_rate
        n_tsave = ceil(Int, T / dt0 / (n_saves - 1))
        dt = T / (n_tsave * (n_saves - 1))
        println("tube S14 $(conf) η=$(eta) α=$(alpha): N0=$(N) dt=$(dt) steps=$(n_tsave * (n_saves - 1))")
        flush(stdout)

        xsaves = zeros(n_saves, Ns_save + 1, 3)
        Qsaves = zeros(n_saves, Ns_save, 3, 3)
        tsaves = zeros(n_saves); Lrest = zeros(n_saves); Nsave = zeros(Int, n_saves)
        Wr_f = zeros(n_saves); Tw = zeros(n_saves); rho_max = zeros(n_saves)
        function record!(nt, t, x, Q, l0s)
            xa = to_array_x(x); Qa = to_array_Q(Q)
            xi, Qi = interpolate_rodstate(xa, Qa, Ns_save)
            xsaves[nt, :, :] = xi; Qsaves[nt, :, :, :] = Qi
            tsaves[nt] = t; Lrest[nt] = sum(l0s); Nsave[nt] = length(Q)
            Wr_f[nt] = writhe_fuller(xa); Tw[nt] = twist_total(Qa)
            rho_max[nt] = maximum(sqrt.(xa[:, 1] .^ 2 .+ xa[:, 2] .^ 2))
        end

        t = 0.0
        record!(1, t, x, Q, l0s)
        for nt = 2:n_saves
            for step = 1:n_tsave
                if step % 5 == 0
                    toobig = l0s .> dxmax
                    if any(toobig)
                        xn, Qn, ln, un, wn = remesh_fast(x, Q, l0s, u0, w, toobig)
                        resize!(x, length(xn)); x .= xn
                        resize!(Q, length(Qn)); Q .= Qn
                        resize!(l0s, length(ln)); l0s .= ln
                        resize!(u0, length(un)); u0 .= un
                        resize!(w, length(wn)); w .= wn
                        rod.N = length(Q)
                        buf = FastBuffers(rod.N)
                    end
                    resize!(cradii, length(x)); fill!(cradii, r_contact)
                    neighbour_points!(nl, arclen, x, cradii)
                    mean_l = sum(l0s) / rod.N
                    rod.gamma_p = gamma_L * mean_l
                    rod.gamma_t = gamma_R * mean_l
                    alpha_node[] = alpha_contact * mean_l
                end
                rk4_step!(x, Q, u0, l0s, dt, rod, buf, extf!, bc!, g_rate, eta)   # growth + remodeling
                t += dt
                sigma == 0 && continue
                for k = 2:rod.N-1               # rotational noise on the tangent
                    amp = sqrt(2 * sigma * dt / l0s[k]) / gamma_R
                    dpsi = amp * (randn(rng) * Q[k][1, :] + randn(rng) * Q[k][2, :])
                    th = norm(dpsi)
                    th == 0.0 && continue
                    u = cross(Q[k][3, :], dpsi) / th
                    Q[k] = Q[k] * rodrigues_s(-th, u)
                end
            end
            record!(nt, t, x, Q, l0s)
        end

        matwrite(out_filename, Dict(
            "xsaves" => xsaves, "Qsaves" => Qsaves, "t" => tsaves, "Lrest" => Lrest, "N" => Nsave,
            "Wr_fuller" => Wr_f, "Tw" => Tw, "rho_max" => rho_max,
            "alpha" => alpha, "eta" => eta, "seed" => seed, "dt" => dt, "confinement" => conf,
            "tube_radius" => tube_radius, "W_wall" => W_wall, "r_contact" => r_contact,
            "alpha_contact" => alpha_contact, "Le" => Le, "g" => g_rate, "B1" => B1, "B3" => B3,
            "sigma" => sigma, "S" => S_stiff))
        return out_filename
    end
end

confs = filter(c -> c in ARGS, ["hard", "free"]); isempty(confs) && (confs = ["hard", "free"])
run_etas = haskey(ENV, "TUBE_ONLY_ETA") ? [parse(Float64, ENV["TUBE_ONLY_ETA"])] : tube_etas
jobs = [(c, eta, ia) for c in confs for eta in run_etas for ia in eachindex(tube_alphas)]
println("---- $(length(jobs)) tube simulations: S14 x $(confs) x η=$(run_etas) x $(length(tube_alphas)) α, to L/Le = $(L_end_tube), σ = $(sigma) ----")
pmap(jobs) do j
    run_tube(j...)
end
println("~success~")
