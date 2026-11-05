## Fast (non-allocating) version of run_simulation_plastic_remeshed (../src/rod_contact2.jl)
#
# Same algorithm, step by step, as the prod driver: RK4 with clamped, time-dependent end
# positions and frames, polyhedral walls, circumradius self-contact, remeshing every 5
# steps (with the per-edge growth weights w_seg), growth and Maxwell remodeling with the
# rest curvature in the material frame. Internal forces and the RK4 sub-steps come from
# ../src/fastkernel.jl. Saves go through prod's interpolate_rodstate, so the output
# has the same layout as prod.
# Checked against run_simulation_plastic_remeshed by validate_fastprod.jl.
#
# Requires ../src/rod_contact2.jl (Plane, interpolate_rodstate) and ../src/fastkernel.jl.

using LinearAlgebra
using StaticArrays

# ---- walls ----------------------------------------------------------------------------
struct FastWalls
    normals::Vector{V3}
    ds::Vector{Float64}
end
FastWalls(walls::Vector{Plane}) = FastWalls([V3(w.normal...) for w in walls], [w.d for w in walls])

function neighbour_walls!(nw::Vector{Vector{Int}}, x::Vector{V3}, W::FastWalls, cutoff::Float64)
    """= make_neighbour_walls: faces with d - n·x <= cutoff, in the wall order"""
    resize!(nw, length(x))
    for i in eachindex(x)
        isassigned(nw, i) ? empty!(nw[i]) : (nw[i] = Int[])
        for f in eachindex(W.ds)
            if W.ds[f] - dot(W.normals[f], x[i]) <= cutoff
                push!(nw[i], f)
            end
        end
    end
    return nw
end

function add_wall_forces!(F::Vector{V3}, x::Vector{V3}, nw::Vector{Vector{Int}}, W::FastWalls, W1::Float64)
    """= ext_forces of the prod scripts: W1 * dist * force_orient on nodes 2..N
    (get_force_components: deepest violated face, scanned in list order)"""
    @inbounds for i = 2:length(x)-1
        min_dist = 0.0
        orient = zero(V3)
        for f in nw[i]
            dist = W.ds[f] - dot(W.normals[f], x[i])
            if dist < min_dist
                min_dist = dist
                orient = -W.normals[f]
            end
        end
        F[i] += W1 * (-min_dist) * (1 + 0 * (-min_dist)) * orient
    end
    return nothing
end

# ---- self-contact ---------------------------------------------------------------------
function neighbour_points!(nl::Vector{Vector{Int}}, arclen::Vector{Float64}, x::Vector{V3},
                           radii::Vector{Float64}; cutoff=1.5, fact=1.0)
    """= make_neighbour_points"""
    n = length(x)
    resize!(nl, n); resize!(arclen, n)
    arclen[1] = 0.0
    for i = 2:n
        arclen[i] = arclen[i-1] + norm(x[i] - x[i-1])
    end
    for i = 1:n
        isassigned(nl, i) ? empty!(nl[i]) : (nl[i] = Int[])
        for j = 1:n
            i == j && continue
            dist = norm(x[i] - x[j])
            if abs(arclen[i] - arclen[j]) > fact * (radii[i] + radii[j])
                if dist <= cutoff * (radii[i] + radii[j])
                    push!(nl[i], j)
                end
            end
        end
    end
    return nl
end

function add_contact_forces!(F::Vector{V3}, x::Vector{V3}, nl::Vector{Vector{Int}},
                             radii::Vector{Float64}, alpha::Float64, fc::Vector{V3})
    """= compute_contact_force_clamped (accumulated in fc, then added to F)"""
    n = length(x)
    resize!(fc, n)
    fill!(fc, zero(V3))
    @inbounds for i = 2:n-1
        isempty(nl[i]) && continue
        dist0 = 1e6
        jc = 0
        for j in nl[i]
            d = norm(x[i] - x[j])
            if d < dist0
                dist0 = d
                jc = j
            end
        end
        overlap = max(radii[i] + radii[jc] - dist0, 0.0)
        orient = normalize(x[jc] - x[i])
        fc[i] += -alpha * overlap * orient
        fc[jc] += alpha * overlap * orient
    end
    fc[1] = zero(V3)
    fc[n] = zero(V3)
    @inbounds for i = 1:n
        F[i] += fc[i]
    end
    return nothing
end

# ---- remeshing ------------------------------------------------------------------------
function remesh_fast(x::Vector{V3}, Q::Vector{M3}, l0s::Vector{Float64}, u0::Vector{V3},
                     w::Vector{Float64}, toobig::BitVector)
    """= remesh (+ split_weights): split the flagged edges in half; the new edge copies the
    parent frame and starts with zero rest curvature at the new interior node"""
    N = length(Q)
    N_new = N + sum(toobig)
    xn = Vector{V3}(undef, N_new + 1)
    Qn = Vector{M3}(undef, N_new)
    ln = zeros(N_new)
    un = zeros(V3, N_new - 1)
    wn = zeros(N_new)
    idx_new = 0
    for idx = 1:N
        idx_new += 1
        if toobig[idx]
            xn[idx_new] = x[idx]
            Qn[idx_new] = Q[idx]
            idx > 1 && (un[idx_new-1] = u0[idx-1])
            ln[idx_new] = 0.5 * l0s[idx]
            wn[idx_new] = 0.5 * w[idx]
            idx_new += 1
            xn[idx_new] = 0.5 * (x[idx] + x[idx+1])
            Qn[idx_new] = Q[idx]
            un[idx_new-1] = zero(V3)
            ln[idx_new] = 0.5 * l0s[idx]
            wn[idx_new] = 0.5 * w[idx]
        else
            xn[idx_new] = x[idx]
            Qn[idx_new] = Q[idx]
            ln[idx_new] = l0s[idx]
            wn[idx_new] = w[idx]
            idx > 1 && (un[idx_new-1] = u0[idx-1])
        end
    end
    xn[end] = x[end]
    return xn, Qn, ln, un, wn
end

# ---- RK4 step with time-dependent bcs and a per-edge growth rate -------------------------
function rk4_step_prod!(t::Float64, x::Vector{V3}, Q::Vector{M3}, u0::Vector{V3}, l0s::Vector{Float64},
                        dt::Float64, rod::FastRod, b::FastBuffers, extf!, bc!,
                        ell_rates::Vector{Float64}, eta::Float64)
    """= timestepRK4plastic_remesh: extf!(F, x) adds walls + contact, bc!(t, x, Q)"""
    F1, F2, F3, F4 = b.F
    T1, T2, T3, T4 = b.T
    x1, x2, x3 = b.xs
    Q1, Q2, Q3 = b.Qs
    forces_torques!(F1, T1, x, Q, u0, l0s, rod, b); extf!(F1, x)
    pos_next!(x1, Q1, x, Q, F1, T1, dt / 2.0, rod); bc!(t, x1, Q1)
    forces_torques!(F2, T2, x1, Q1, u0, l0s, rod, b); extf!(F2, x1)
    pos_next!(x2, Q2, x, Q, F2, T2, dt / 2.0, rod); bc!(t + dt / 2.0, x2, Q2)
    forces_torques!(F3, T3, x2, Q2, u0, l0s, rod, b); extf!(F3, x2)
    pos_next!(x3, Q3, x, Q, F3, T3, dt, rod); bc!(t + dt / 2.0, x3, Q3)
    forces_torques!(F4, T4, x3, Q3, u0, l0s, rod, b); extf!(F4, x3)
    @inbounds for i in eachindex(F1)
        F1[i] = F1[i] + 2.0 * F2[i] + 2.0 * F3[i] + F4[i]
    end
    @inbounds for i in eachindex(T1)
        T1[i] = T1[i] + 2.0 * T2[i] + 2.0 * T3[i] + T4[i]
    end
    pos_next!(x, Q, x, Q, F1, T1, dt / 6.0, rod); bc!(t + dt, x, Q)
    @inbounds for i in eachindex(l0s)
        l0s[i] = l0s[i] + ell_rates[i] * dt
    end
    @inbounds for i = 1:rod.N-1
        voro = (l0s[i] + l0s[i+1]) / 2.0
        kap = rotlog2_s(Q[i]' * Q[i+1]) / voro
        target = rod.material_u0 ? Q[i] * kap : kap
        u0[i] = (u0[i] + eta * dt * target) / (1 + eta * dt)
    end
    return nothing
end

function fractional_arclength!(sc::Vector{Float64}, x::Vector{V3})
    """= get_fractional_arclength"""
    n = length(x)
    resize!(sc, n)
    sc[1] = 0.0
    for i = 2:n
        sc[i] = sc[i-1] + norm(x[i] - x[i-1])
    end
    tot = sc[n]
    for i = 1:n
        sc[i] /= tot
    end
    return sc
end

# ---- driver: = run_simulation_plastic_remeshed (proba = 0, elastic = false) -----------
function run_fastprod(xpos0::Matrix{Float64}, Qs0::Array{Float64,3}, dt::Float64, Ns::Int, T::Float64,
                      n_saves::Int, S::Float64, St::Float64, B::Float64, Bt::Float64,
                      gamma_p::Float64, gamma_t::Float64, bc!, circum_radii, alpha::Float64,
                      walls::Vector{Plane}, cutoff::Float64, W1::Float64, growth_rates,
                      plastic_rate::Float64, dxmax::Float64; material_u0::Bool=true, show_progress::Bool=false)
    """bc!(t, x, Q) imposes the clamped ends on Vector{V3} / Vector{M3} in place"""
    N = size(Qs0, 1)
    xsaves = zeros(n_saves, Ns + 1, 3)
    Qsaves = zeros(n_saves, Ns, 3, 3)
    n_tsave = round(Int, T / dt / n_saves)

    x_i, Q_i = interpolate_rodstate(xpos0, Qs0, Ns)
    xsaves[1, :, :] = x_i
    Qsaves[1, :, :, :] = Q_i
    t = 0.0

    x = to_static_x(xpos0)
    Q = [rownormalize_s(M3(Qs0[k, :, :])) for k = 1:N]
    l0s = [norm(x[i+1] - x[i]) for i = 1:N]
    u0 = zeros(V3, N - 1)
    for i = 1:N-1 # = get_u0 (material frame)
        voro = (l0s[i] + l0s[i+1]) / 2.0
        kap = rotlog2_s(Q[i]' * Q[i+1]) / voro
        u0[i] = material_u0 ? Q[i] * kap : kap
    end
    w_seg = ones(N)

    rod = FastRod(N, S, St, B, Bt, gamma_p, gamma_t; material_u0=material_u0)
    buf = FastBuffers(N)
    W = FastWalls(walls)
    nw = Vector{Vector{Int}}(); nl = Vector{Vector{Int}}()
    arclen = Float64[]; sc = Float64[]; cradii = Float64[]; ell_rates = Float64[]; fc = V3[]

    function refresh_rates!(t)
        fractional_arclength!(sc, x)
        resize!(cradii, length(x)); resize!(ell_rates, length(x) - 1)
        for i in eachindex(sc)
            cradii[i] = circum_radii(t, sc[i])
        end
        for i = 1:length(x)-1 # growth sampled at the end node of each edge, as in prod
            ell_rates[i] = w_seg[i] * growth_rates(t, sc[i+1])
        end
    end
    function extf!(F, xx)
        add_wall_forces!(F, xx, nw, W, W1)
        add_contact_forces!(F, xx, nl, cradii, alpha, fc)
    end

    # first time point
    neighbour_walls!(nw, x, W, cutoff)
    refresh_rates!(t)
    neighbour_points!(nl, arclen, x, cradii)
    rk4_step_prod!(t, x, Q, u0, l0s, dt, rod, buf, extf!, bc!, ell_rates, plastic_rate)

    for nt = 1:(n_saves-1)
        for step = 1:n_tsave
            t += dt
            if step % 5 == 0
                toobig = l0s .> dxmax
                if any(toobig)
                    # update in place: x and w_seg are captured by the closures above
                    xn, Qn, ln, un, wn = remesh_fast(x, Q, l0s, u0, w_seg, toobig)
                    resize!(x, length(xn)); x .= xn
                    resize!(Q, length(Qn)); Q .= Qn
                    resize!(l0s, length(ln)); l0s .= ln
                    resize!(u0, length(un)); u0 .= un
                    resize!(w_seg, length(wn)); w_seg .= wn
                    rod.N = length(Q)
                    buf = FastBuffers(rod.N)
                end
                neighbour_walls!(nw, x, W, cutoff)
                refresh_rates!(t)
                neighbour_points!(nl, arclen, x, cradii)
            else
                refresh_rates!(t)
            end
            rk4_step_prod!(t, x, Q, u0, l0s, dt, rod, buf, extf!, bc!, ell_rates, plastic_rate)
        end
        x_i, Q_i = interpolate_rodstate(to_array_x(x), to_array_Q(Q), Ns)
        xsaves[nt+1, :, :] = x_i
        Qsaves[nt+1, :, :, :] = Q_i
        show_progress && (print("\r save $(nt+1)/$(n_saves)  N=$(rod.N)"); flush(stdout))
    end
    return xsaves, Qsaves
end
