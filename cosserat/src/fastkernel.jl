## Non-allocating copy of the prod rod kernel (../src/rod_contact2.jl) for the reduced/ runs
#
# Same arithmetic, term by term, as
#   internal_forces_and_torques(x, Q, rod, u0, l0s), pos_next, timestepRK4plastic_remesh
#   and get_u0_plastic
# specialized to: uniform moduli, no walls, no self-contact, no external torques, and an
# external force callback. Positions are Vector{SVector{3}} and frames Vector{SMatrix{3,3}}
# (rows = directors, as in prod). Checked against the prod kernel by validate_fast.jl.
#
# material_u0 = true (default) stores u0 in the material frame, Qᵢ κ, so that the elastic
#   energy is invariant under rigid rotations also when u0 ≠ 0 (remodeling runs). This is
#   the convention of ../src/rod_contact2.jl since the 2026-10-19 fix.
# material_u0 = false reproduces the pre-fix kernel (../src/rod_contact2_labu0.jl), where u0
#   is stored in the same lab-frame form as rotLog_2(Qᵢ'Qᵢ₊₁)/ℓ̄ and does not rotate with
#   the material.

using StaticArrays
using LinearAlgebra

const V3 = SVector{3,Float64}
const M3 = SMatrix{3,3,Float64,9}
const EZ3 = V3(0.0, 0.0, 1.0)

skew_s(u) = M3(0.0, u[3], -u[2], -u[3], 0.0, u[1], u[2], -u[1], 0.0)   # = skewVec2Mat(u)
unskew_s(U) = V3(U[3, 2], -U[3, 1], U[2, 1])                            # = skewMat2Vec(U)

function rodrigues_s(theta::Float64, u::V3)
    U = skew_s(u)
    return one(M3) + sin(theta) * U + (1 - cos(theta)) * U * U
end

function rotlog2_s(M::M3)
    """= rotLog_2(M)"""
    u = 0.5 * unskew_s(M - M')
    n = norm(u)
    n == 0.0 && return zero(V3)
    u = u / n
    theta = atan(-tr(skew_s(u) * M) / 2.0, (tr(M) - 1.0) / 2.0)
    return theta * u
end

function rownormalize_s(M::M3)
    d = V3(norm(M[1, :]), norm(M[2, :]), norm(M[3, :]))
    return M ./ d
end

mutable struct FastRod
    N::Int
    Smat::M3
    Bmat::M3
    gamma_p::Float64
    gamma_t::Float64
    material_u0::Bool
    twist_drag::Float64   # rotational drag about d3 = twist_drag * gamma_t (1: isotropic, as in prod)
end

FastRod(N, S, St, B, Bt, gamma_p, gamma_t; material_u0=true, twist_drag=1.0) =
    FastRod(N, M3(S, 0, 0, 0, S, 0, 0, 0, St), M3(B, 0, 0, 0, B, 0, 0, 0, Bt), gamma_p, gamma_t,
            material_u0, twist_drag)

struct FastBuffers
    len::Vector{Float64}
    ten::Vector{V3}
    mom::Vector{V3}
    mom2::Vector{V3}
    F::NTuple{4,Vector{V3}}
    T::NTuple{4,Vector{V3}}
    xs::NTuple{3,Vector{V3}}
    Qs::NTuple{3,Vector{M3}}
end

function FastBuffers(N::Int)
    vN1() = zeros(V3, N + 1)
    vN() = zeros(V3, N)
    qN() = fill(one(M3), N)
    return FastBuffers(zeros(N), vN(), zeros(V3, N - 1), zeros(V3, N - 1),
                       (vN1(), vN1(), vN1(), vN1()), (vN(), vN(), vN(), vN()),
                       (vN1(), vN1(), vN1()), (qN(), qN(), qN()))
end

# conversions to / from the prod array layout
to_static_x(x::AbstractMatrix) = [V3(x[i, 1], x[i, 2], x[i, 3]) for i = 1:size(x, 1)]
to_static_Q(Q::AbstractArray) = [M3(Q[k, :, :]) for k = 1:size(Q, 1)]
to_static_u(u::AbstractMatrix) = [V3(u[i, 1], u[i, 2], u[i, 3]) for i = 1:size(u, 1)]
to_array_x(x::Vector{V3}) = permutedims(reduce(hcat, x))
function to_array_Q(Q::Vector{M3})
    A = zeros(length(Q), 3, 3)
    for k in eachindex(Q)
        A[k, :, :] = Q[k]
    end
    return A
end

function forces_torques!(F::Vector{V3}, Tq::Vector{V3}, x::Vector{V3}, Q::Vector{M3},
                         u0::Vector{V3}, l0s::Vector{Float64}, rod::FastRod, b::FastBuffers)
    """= internal_forces_and_torques(x, Q, rod, u0, l0s)"""
    N = rod.N
    @inbounds for i = 1:N
        b.len[i] = norm(x[i+1] - x[i])
    end
    fill!(Tq, zero(V3))
    @inbounds for i = 1:N-1 # interior vertices
        ri = x[i+1] - x[i]
        li = b.len[i]
        voro0 = (l0s[i] + l0s[i+1]) / 2.0
        voroD = (li + b.len[i+1]) / 2.0
        Eps = voroD / voro0
        ti = ri / li
        elo = li / l0s[i]
        qs = Q[i]
        sig = rod.Smat * (qs * elo * ti - EZ3)
        kap = rotlog2_s(qs' * Q[i+1]) / voro0
        b.ten[i] = qs' * sig
        if rod.material_u0
            m_mat = rod.Bmat * (qs * kap - u0[i]) / Eps^3
        else
            m_mat = rod.Bmat * (qs * (kap - u0[i])) / Eps^3
        end
        b.mom[i] = qs' * m_mat
        b.mom2[i] = qs' * cross(qs * kap, m_mat) * voro0
        Tq[i] += cross(qs' * sig, ti) * l0s[i]
    end
    # last edge
    ri = x[N+1] - x[N]
    li = b.len[N]
    ti = ri / li
    elo = li / l0s[N]
    sig = rod.Smat * (Q[N] * elo * ti - EZ3)
    b.ten[N] = Q[N]' * sig
    Tq[N] += cross(Q[N]' * sig, ti) * l0s[N]
    # forces = delta_op(tensions, N)
    F[1] = b.ten[1]
    @inbounds for j = 2:N
        F[j] = b.ten[j] - b.ten[j-1]
    end
    F[N+1] = -b.ten[N]
    # torques += delta_op(moments, N-1) + A_op(moments2, N-1)
    Tq[1] += b.mom[1] + 0.5 * b.mom2[1]
    @inbounds for j = 2:N-1
        Tq[j] += (b.mom[j] - b.mom[j-1]) + 0.5 * (b.mom2[j] + b.mom2[j-1])
    end
    Tq[N] += -b.mom[N-1] + 0.5 * b.mom2[N-1]
    return nothing
end

function pos_next!(xo::Vector{V3}, Qo::Vector{M3}, x::Vector{V3}, Q::Vector{M3},
                   F::Vector{V3}, Tq::Vector{V3}, dt::Float64, rod::FastRod)
    """= pos_next (xo, Qo may alias x, Q)"""
    @inbounds for i in eachindex(x)
        xo[i] = x[i] + dt * F[i] / rod.gamma_p
    end
    @inbounds for i = 1:rod.N
        angvel = Tq[i] / rod.gamma_t
        if rod.twist_drag != 1.0 # anisotropic drag: rotation about d3 has drag twist_drag * gamma_t
            d3 = Q[i][3, :]
            angvel += (1 / rod.twist_drag - 1) * dot(angvel, d3) * d3
        end
        theta = norm(angvel)
        Qn = theta == 0.0 ? Q[i] : Q[i] * rodrigues_s(theta * dt, angvel / theta)
        Qo[i] = rownormalize_s(Qn)
    end
    return nothing
end

function rk4_step!(x::Vector{V3}, Q::Vector{M3}, u0::Vector{V3}, l0s::Vector{Float64},
                   dt::Float64, rod::FastRod, b::FastBuffers, extf!, bc!, g_rate::Float64, eta::Float64)
    """= timestepRK4plastic_remesh with ell_rates = g_rate .* l0s and plastic_rate = eta,
    updating (x, Q, u0, l0s) in place. extf!(F, x) adds the external forces to F,
    bc!(x, Q) imposes the boundary conditions."""
    F1, F2, F3, F4 = b.F
    T1, T2, T3, T4 = b.T
    x1, x2, x3 = b.xs
    Q1, Q2, Q3 = b.Qs
    # stage 1
    forces_torques!(F1, T1, x, Q, u0, l0s, rod, b); extf!(F1, x)
    pos_next!(x1, Q1, x, Q, F1, T1, dt / 2.0, rod); bc!(x1, Q1)
    # stage 2
    forces_torques!(F2, T2, x1, Q1, u0, l0s, rod, b); extf!(F2, x1)
    pos_next!(x2, Q2, x, Q, F2, T2, dt / 2.0, rod); bc!(x2, Q2)
    # stage 3
    forces_torques!(F3, T3, x2, Q2, u0, l0s, rod, b); extf!(F3, x2)
    pos_next!(x3, Q3, x, Q, F3, T3, dt, rod); bc!(x3, Q3)
    # stage 4
    forces_torques!(F4, T4, x3, Q3, u0, l0s, rod, b); extf!(F4, x3)
    # result
    @inbounds for i in eachindex(F1)
        F1[i] = F1[i] + 2.0 * F2[i] + 2.0 * F3[i] + F4[i]
    end
    @inbounds for i in eachindex(T1)
        T1[i] = T1[i] + 2.0 * T2[i] + 2.0 * T3[i] + T4[i]
    end
    pos_next!(x, Q, x, Q, F1, T1, dt / 6.0, rod); bc!(x, Q)
    # growth, then plasticity at the new rest lengths (= get_u0_plastic)
    @inbounds for i in eachindex(l0s)
        l0s[i] = l0s[i] + g_rate * l0s[i] * dt
    end
    @inbounds for i = 1:rod.N-1
        voro = (l0s[i] + l0s[i+1]) / 2.0
        kap = rotlog2_s(Q[i]' * Q[i+1]) / voro
        target = rod.material_u0 ? Q[i] * kap : kap
        u0[i] = (u0[i] + eta * dt * target) / (1 + eta * dt)
    end
    return nothing
end
