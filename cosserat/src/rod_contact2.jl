#using DifferentialEquations
using LinearAlgebra
#using GeometryBasics
#using ComponentArrays
using ProgressMeter
using Interpolations

#import Base.:*
#import Base.:+


include("./walls.jl")
include("./contact.jl")

struct RodElements
    rest_l::Vector{Float64}
    voro_0::Vector{Float64}
    gamma_p::Float64
    gamma_t::Float64
    shearMat::AbstractArray{Float64}
    moduli::AbstractArray{Float64}
    B::Vector{Float64}
    Bt::Vector{Float64}
    u0::AbstractArray{Float64}
    N::Int
end



function get_torsion_2(d3_vec)
    N = length(d3_vec)
    tw = zeros(N)
    tw[1] = 0.
    for i=2:N-2
        t1 = d3_vec[i-1]
        t2 = d3_vec[i]
        t3 = d3_vec[i+1]
        axis1 = cross(t1, t2)
        axis1 = axis1/norm(axis1)
        axis2 = cross(t2, t3)
        axis2 = axis2/norm(axis2)
        oAxis = cross(t2, axis1)
        cosT = dot(axis2, axis1)
        sinT = dot(axis2, oAxis)
        tw[i] = atan(sinT, cosT)
    end
    tw[end-1] = tw[end-2]
    tw[end] = tw[end-1]
    return tw
end


function get_twistdensity(d1_vec, d2_vec, d3_vec)
    N = length(d1_vec)
    tw = zeros(N)
    tw[1] = 0.
    u1 = zeros(N,3)
    u2 = zeros(N,3)
    u1[1,:] = d1_vec[1]
    u2[2,:] = d2_vec[1]
    for i=2:N-1
        t1 = d3_vec[i-1]
        t2 = d3_vec[i]
        axis = cross(t1, t2)
        if norm(axis) < 1e-4
            axis =  d1_vec[i-1]
        else
            axis = axis/norm(axis)
        end
        oAxis = cross(t2, axis)
        cosA = dot(t1, t2)
        sinA = dot(t1, oAxis)
        #tw[i] = atan(sinA, cosA)
        u1[i, :] = cosA * u1[i-1, :] + cross(axis, u1[i-1,:]) * sinA + axis * dot(axis, u1[i-1,:]) * (1.0 -cosA)
        u2[i, :] = cosA * u2[i-1, :] + cross(axis, u2[i-1,:]) * sinA + axis * dot(axis, u2[i-1,:]) * (1.0 -cosA)
        u1[i,:] /= norm(u1[i,:])
        u2[i,:] /= norm(u2[i,:])
        tw[i] = atan(dot(d1_vec[i], u2[i,:]), dot(d1_vec[i], u1[i,:]))
    end
    tw[end] = tw[end-1]
    return tw
end

function rotation_angle(v1::Vector{Float64}, v2::Vector{Float64})
	# BUGFIX: parentheses — take acos of the normalized dot product, not the other way around
	return acos(dot(v1, v2)/(norm(v1)*norm(v2)))
end

function skewMat2Vec(U::AbstractArray{Float64})
    return [U[3,2], - U[3,1], U[2,1]]
end

function skewVec2Mat(u::Vector{Float64})
    return [ 0.0 -u[3]  u[2];
             u[3]  0.0 -u[1];
            -u[2] u[1]   0.0]
end

function Rodrigues(theta::Float64, u::Vector{Float64})
    U = skewVec2Mat(u)
    return I + sin(theta) * U + (1-cos(theta)) * U*U
end

function RodriguesMat(M::Matrix{Float64},theta::Float64, u::Vector{Float64} )
    return M * Rodrigues(theta, u)
end

function rotLog(Mat::AbstractArray{Float64})
    """Returns the matrix logarithm of a rotation matrix"""
    theta = acos((tr(Mat)-1.0)/2.0)
    if theta == 0
        return [0., 0., 0.]
    else
        return .5*theta/sin(theta) * skewMat2Vec(Mat - Mat')
    end
end

function rotLog_2(Mat::AbstractArray{Float64})
    """Returns the matrix logarithm of a rotation matrix"""
    #theta = acos((tr(Mat)-1.0)/2.0)
    u = 0.5*skewMat2Vec(Mat - Mat')
    if norm(u) == 0.
        return [0., 0., 0.]
    else
        u = u / norm(u)
        theta = atan(-tr(skewVec2Mat(u)*Mat)/2., (tr(Mat)-1.0)/2.)
        return theta * u
    end
end

function rotVec(A, B)
    """finds the rotation vector v that rotates rotation matrix A into B
    e.g. R(v)A = B
    """
    d31 = A[3,:]/norm(A[3,:])
    d32 = B[3,:]/norm(B[3,:])
    k = cross(d31, d32)
    kp = cross(k, d31)
    phi = atan(dot(kp, d32), dot(d31, d32))
    return phi * k
end

function delta_op(x::AbstractArray{Float64}, N::Int)
    Dx = zeros(N+1, 3)
    Dx[1,:] = x[1,:]
    Dx[2:end-1,:] = x[2:end,:] - x[1:end-1,:]
    Dx[end,:] = -x[end,:]
    return Dx
end

function A_op(x::AbstractArray{Float64}, N::Int)
    Ax = zeros(N+1, 3)
    Ax[1,:] = 0.5*x[1,:]
    Ax[2:end-1,:] = 0.5*(x[2:end,:] + x[1:end-1,:])
    Ax[end,:] = 0.5*x[end,:]
    return Ax
end

function get_adaptedQ(xpos::AbstractArray{Float64}, Tw::Real)
    N = size(xpos, 1)-1
    Qs = zeros(N,3,3)
    # permute vector to get a candidate for orthogonal vector to d3
    r3 = xpos[2,:] - xpos[1,:]
    d3 = r3 / norm(r3)
    u1 = [r3[3], r3[1], r3[2]]
    u1 = u1 - dot(u1, d3)*d3
    d1 = u1 / norm(u1)
    d2 = cross(d3, d1)
    Qs[1,1,:] = d1
    Qs[1,2,:] = d2
    Qs[1,3,:] = d3
    for i=2:N
        r3 = xpos[i+1,:] - xpos[i,:]
        d3 = r3 / norm(r3)
        # find matrix rotating Qs_i into Qs_i+1
        u = cross(Qs[i-1,3,:], d3)
        if norm(u) < 1e-12
            Qs[i,:,:] = Qs[i-1,:,:]
        else
            u = u/norm(u)
            ux = cross(Qs[i-1,3,:], u)
            theta = atan(dot(d3 ,ux), dot(d3, Qs[i-1,3,:]))
            R = Rodrigues(theta, u)
            Qs[i,:,:] .= Qs[i-1,:,:] * R
        end
        # add in twist : rotate Qs_i around d3 by Tw/N
        R_tw = Rodrigues(Tw/(N-1), Qs[i,3,:])
        Qs[i,:,:] .= Qs[i,:,:] * R_tw
    end
    return Qs
end


function make_Rod(xpos::AbstractArray{Float64}, S::Float64, St::Float64, B::Float64, Bt::Float64, gamma_p::Float64, gamma_t::Float64)
    """returns rest legnths, rest voronoi domains, rest curvatures, stiffness and shear matrices for a rod of shear rigidity S, stretch modulus St, bending modulus B, twist modulus Bt"""

    N = size(xpos,1)-1

    moduli = zeros(N,3,3)
    shearMat = zeros(N,3,3)
    u0 = zeros(N,3)
    nu = 0.0

    rest_l = zeros(N)
    voro_0 = zeros(N-1)

    #I1 = pi/4 *r^4  #1.0 /12.0 * (3.0 * r^2 + rest_l[i]^2)
    #I3 = 2*I1

    for i = 1:N
        rest_l[i] = norm(xpos[i,:]-xpos[i+1,:])
        shearMat[i,:,:] = diagm([S, S, St])
        moduli[i,:,:] = diagm([B, B, Bt])
    end
    voro_0 = (rest_l[1:(N-1)] + rest_l[2:N]) / 2.0

    B_vec = fill(B, N)
    Bt_vec = fill(Bt, N)

    rod = RodElements(rest_l,voro_0,
                        gamma_p,gamma_t,
                        shearMat,moduli,
                        B_vec, Bt_vec,
                        u0, N)
    return rod
end

function make_Rod(xpos::AbstractArray{Float64}, S::Float64, St::Float64, B::Vector{Float64}, Bt::Vector{Float64}, gamma_p::Float64, gamma_t::Float64)
    """returns rest legnths, rest voronoi domains, rest curvatures, stiffness and shear matrices for a rod of shear rigidity S, stretch modulus St, bending modulus B, twist modulus Bt"""

    N = size(xpos,1)-1
    moduli = zeros(N,3,3)
    shearMat = zeros(N,3,3)
    u0 = zeros(N,3)

    rest_l = zeros(N)
    voro_0 = zeros(N-1)

    B_vec = zeros(N)
    Bt_vec = zeros(N)

    for i = 1:N
        rest_l[i] = norm(xpos[i,:]-xpos[i+1,:])
        shearMat[i,:,:] = diagm([S, S, St])
        moduli[i,:,:] = diagm([B[i], B[i], Bt[i]])
        B_vec[i] = B[i]
        Bt_vec[i] = Bt[i]
    end
    voro_0 = (rest_l[1:(N-1)] + rest_l[2:N]) / 2.0

    rod = RodElements(rest_l,voro_0,
                        gamma_p,gamma_t,
                        shearMat,moduli,
                        B_vec, Bt_vec,
                        u0, N)
    return rod
end

function make_Rod(xpos::AbstractArray{Float64}, S::Vector{Float64}, St::Vector{Float64}, B::Vector{Float64}, Bt::Vector{Float64}, gamma_p::Float64, gamma_t::Float64)
    """returns rest legnths, rest voronoi domains, rest curvatures, stiffness and shear matrices for a rod of shear rigidity S, stretch modulus St, bending modulus B, twist modulus Bt"""

    N = size(xpos,1)-1
    moduli = zeros(N,3,3)
    shearMat = zeros(N,3,3)
    u0 = zeros(N,3)

    rest_l = zeros(N)
    voro_0 = zeros(N-1)

    B_vec = zeros(N+1)
    Bt_vec = zeros(N+1)

    for i = 1:N
        rest_l[i] = norm(xpos[i,:]-xpos[i+1,:])
        shearMat[i,:,:] = diagm([S[i], S[i], St[i]])
        moduli[i,:,:] = diagm([B[i], B[i], Bt[i]])
        B_vec[i] = B[i]
        Bt_vec[i] = Bt[i]
    end
    voro_0 = (rest_l[1:(N-1)] + rest_l[2:N]) / 2.0

    rod = RodElements(rest_l,voro_0,
                        gamma_p,gamma_t,
                        shearMat,moduli,
                        B_vec, Bt_vec,
                        u0, N)
    return rod
end


function set_curvature(rod, u0_new)
    """Returns a new RodElements with updated rest curvature u0_new"""
    # BUGFIX: RodElements has no E/nu fields; carry over the actual fields
    # (shearMat, moduli, B, Bt) instead.
    rod_new = RodElements(rod.rest_l, rod.voro_0,
                        rod.gamma_p, rod.gamma_t,
                        rod.shearMat, rod.moduli,
                        rod.B, rod.Bt,
                        u0_new, rod.N)
    return rod_new
end

function set_curvature!(rod::RodElements, u0::AbstractArray{Float64})
    """In-place update of rest curvature u0"""
    rod.u0 .= u0
end

function get_fractional_arclength(x)
    """Returns the fractional arclengths along the rod defined by positions x"""
    N = size(x,1)-1
    arclengths = zeros(N)
    for i = 1:N
        arclengths[i] = norm(x[i+1,:] - x[i,:])
    end
    total_length = sum(arclengths)
    fractional_arclengths = [0.; cumsum(arclengths)] / total_length
    return fractional_arclengths
end


# gamma *v = internal_forces + external_forces
# gamma * omega  = internal_torques + external_forces
# d_t x = v
# d_t Q = omega x Q # represent that for dt via exact integration for constant omega/ matrix exp

function internal_forces_and_torques(xpos::AbstractArray{Float64}, Qs::AbstractArray{Float64},
                             rod::RodElements)

    #
    tensions = zeros(rod.N,3)
    moments = zeros(rod.N-1,3)
    moments2 = zeros(rod.N-1,3)
    forces = zeros(rod.N+1,3)
    torques = zeros(rod.N,3)

    lengths_l = zeros(rod.N)
    # compute dilatations
    @inbounds for i=1:rod.N
        ri = xpos[i+1,:]-xpos[i,:]
        lengths_l[i] = norm(ri)
    end

    @inbounds for i = 1:(rod.N-1) # interior vertices
        ri = xpos[i+1,:]-xpos[i,:]
        li = lengths_l[i]
        voro_D = (li + lengths_l[i+1])/2.0
        Eps = voro_D/ rod.voro_0[i]
        ti = ri / li
        elo = li / rod.rest_l[i]
        # trick: use property that Q*d3 = [0., 0., 1.] constant
        qs = @view Qs[i,:,:]
        Shear_dot_sigma_i = rod.shearMat[i,:,:] * (qs * elo*ti - [0., 0., 1.0])
        kappa_i = rotLog_2(qs' * Qs[i+1,:,:]) / rod.voro_0[i]
        tensions[i,:] = qs' *  Shear_dot_sigma_i
        # BUGFIX: the stiffness B = diag(B,B,Bt) is diagonal only in the material
        # frame, but kappa_i is a lab-frame curvature vector. Rotate the curvature
        # deviation into the material frame (qs*...), apply the moduli there, then
        # rotate the resulting couple back to the lab frame (qs'*...). This reduces
        # to the previous code exactly when B is isotropic (Bt == B), but for
        # Bt != B it correctly applies the twist modulus about the local d3 rather
        # than the fixed lab z-axis.
        # BUGFIX (2026-10-19): the rest curvature u0 is stored in the MATERIAL frame
        # (see get_u0), so only kappa_i is rotated. Previously u0 was stored in the same
        # lab-frame form as kappa_i and did not rotate with the rod, which broke rotation
        # invariance of the elastic energy whenever u0 != 0.
        m_mat = rod.moduli[i,:,:] * (qs * kappa_i - rod.u0[i,:]) / Eps^3
        moments[i,:] = qs' * m_mat
        moments2[i,:] = qs' * cross(qs * kappa_i, m_mat) * rod.voro_0[i]
        torques[i,:] += cross(qs' *  Shear_dot_sigma_i, ti) *rod.rest_l[i]
    end
    # last edge
    ri = xpos[rod.N+1,:]-xpos[rod.N,:]
    li = lengths_l[rod.N]
    ti = ri / li
    elo = li / rod.rest_l[rod.N]
    Shear_dot_sigma_i = rod.shearMat[rod.N,:,:] * (Qs[rod.N,:,:] * elo*ti - [0., 0., 1.0])
    tensions[rod.N,:] = Qs[rod.N,:,:]' * Shear_dot_sigma_i
    torques[rod.N,:] +=  cross(Qs[rod.N,:,:]' * Shear_dot_sigma_i, ti) *rod.rest_l[rod.N]
    forces = delta_op(tensions, rod.N) # force is (N+1)x3

    torques += delta_op(moments, rod.N-1) + A_op(moments2, rod.N-1) # torque is Nx3

    return forces, torques

end


function internal_forces_and_torques(xpos::AbstractArray{Float64}, Qs::AbstractArray{Float64},
    rod::RodElements, u0::AbstractArray{Float64}, l0s::AbstractArray{Float64})

    #
    tensions = zeros(rod.N,3)
    moments = zeros(rod.N-1,3)
    moments2 = zeros(rod.N-1,3)
    forces = zeros(rod.N+1,3)
    torques = zeros(rod.N,3)

    lengths_l = zeros(rod.N)
    voro_0 = zeros(rod.N-1)
    # compute dilatations
    @inbounds for i=1:rod.N
        ri = xpos[i+1,:]-xpos[i,:]
        lengths_l[i] = norm(ri)
    end

    @inbounds for i = 1:(rod.N-1) # interior vertices
        ri = xpos[i+1,:]-xpos[i,:]
        li = lengths_l[i]
        voro_0[i] = (l0s[i]+l0s[i+1])/2.0
        voro_D = (li + lengths_l[i+1])/2.0
        Eps = voro_D/ voro_0[i]
        ti = ri / li
        elo = li / l0s[i]
        # trick: use property that Q*d3 = [0., 0., 1.] constant
        qs = @view Qs[i,:,:]
        Shear_dot_sigma_i = rod.shearMat[i,:,:] * (qs * elo*ti - [0., 0., 1.0])
        #println(qs * Shear_dot_sigma_i)
        kappa_i = rotLog_2(qs' * Qs[i+1,:,:]) / voro_0[i]
        #kappa_i = rotVec(Qs[i+1,:,:], qs) / rod.voro_0[i]
        #println(kappa_i)
        tensions[i,:] = qs' *  Shear_dot_sigma_i
        # BUGFIX: apply the moduli in the material frame and rotate the couple back
        # to the lab frame (see the matching comment in the other method). kappa_i is
        # a lab-frame curvature but B = diag(B,B,Bt) is diagonal only in the material
        # frame; equivalent to the old code when B is isotropic.
        # BUGFIX (2026-10-19): u0 is stored in the material frame (see get_u0 and the
        # matching comment in the other method), so only kappa_i is rotated.
        m_mat = rod.moduli[i,:,:] * (qs * kappa_i - u0[i,:]) / Eps^3
        moments[i,:] = qs' * m_mat
        moments2[i,:] = qs' * cross(qs * kappa_i, m_mat) * voro_0[i]
        torques[i,:] += cross(qs' *  Shear_dot_sigma_i, ti) *l0s[i]
    end
    # last edge
    ri = xpos[rod.N+1,:]-xpos[rod.N,:]
    li = lengths_l[rod.N]
    ti = ri / li
    elo = li / l0s[rod.N]
    #println(elo)
    Shear_dot_sigma_i = rod.shearMat[rod.N,:,:] * (Qs[rod.N,:,:] * elo*ti - [0., 0., 1.0])
    tensions[rod.N,:] = Qs[rod.N,:,:]' * Shear_dot_sigma_i
    torques[rod.N,:] +=  cross(Qs[rod.N,:,:]' * Shear_dot_sigma_i, ti) * l0s[rod.N]
    #println(Qs[N,:,:]'*Shear_dot_sigma_i)
    forces = delta_op(tensions, rod.N) # force is (N+1)x3

    torques += delta_op(moments, rod.N-1) + A_op(moments2, rod.N-1) # torque is Nx3

    return forces, torques

end

function get_u0(xpos::Matrix{Float64}, Qs::Array{Float64, 3}, l0s::Vector{Float64})
    """Return u0 adapted to the current x, Q, l0s, in the material frame of edge i"""
    N = size(xpos,1)-1
    u0 = zeros(N-1,3)
    for i =1:(N-1)
        voro = (l0s[i]+l0s[i+1])/2.0
        kappa_i = rotLog_2(Qs[i,:,:]' * Qs[i+1,:,:]) / voro
        u0[i,:] = Qs[i,:,:] * kappa_i   # BUGFIX (2026-10-19): material frame
    end
    return u0
end

function get_l0s(xpos::Matrix{Float64}, Qs::Array{Float64, 3})
    """Return l0s adapted to the current x, Q"""
    N = size(xpos,1)-1
    l0s = zeros(N)
    for i =1:N
        l0s[i] = norm(xpos[i+1,:]-xpos[i,:])
    end
    return l0s
end


function get_u0_plastic(xpos::Matrix{Float64}, Qs::Array{Float64, 3}, l0s::Vector{Float64}, u0::Matrix{Float64}, dt::Float64, coeff::Float64)
    """Maxwell relaxation of u0 towards the current curvature, in the material frame"""
    N = size(xpos,1)-1
    u0n = zeros(N-1,3)
    for i =1:(N-1)
        voro = (l0s[i]+l0s[i+1])/2.0
        kappa_i = rotLog_2(Qs[i,:,:]' * Qs[i+1,:,:]) / voro
        # BUGFIX (2026-10-19): relax towards the material-frame curvature Q_i*kappa_i
        u0n[i,:] =  (u0[i,:] + coeff*dt*(Qs[i,:,:] * kappa_i))/(1+coeff*dt)
    end
    return u0n
end


function pos_next(xpos::Matrix{Float64}, Qs::Array{Float64, 3}, forces::AbstractArray{Float64}, torques::AbstractArray{Float64}, dt::Float64, gamma_p::Float64, gamma_t::Float64, N::Int)
    xnext = xpos + dt*forces/gamma_p # new positions
    Qnext = zeros(N, 3,3)
    for i=1:N
        angvel = torques[i,:]/gamma_t
        theta = norm(angvel)
        if theta == 0.
            Qnext[i,:,:] = Qs[i,:,:]
        else
            Qnext[i,:,:] = RodriguesMat(Qs[i,:,:], theta*dt, angvel/theta)
        end
        Qnext[i,1,:] .= Qnext[i,1,:] / norm(Qnext[i,1,:])
        Qnext[i,2,:] .= Qnext[i,2,:] / norm(Qnext[i,2,:])
        Qnext[i,3,:] .= Qnext[i,3,:] / norm(Qnext[i,3,:])
    end
    return xnext, Qnext
end

function free_bcs()
    """ Generates a boundary condition function for free boundary conditions"""
    function f1(t, x, Q)
        return x, Q
    end
    return f1
end

function free_ext_forces()
    """ Generates a trivial external force function """
    function f1(t, x, Q)
        # BUGFIX: zeros(::Matrix) is invalid; return a zero array of matching shape
        return zero(x)
    end
    return f1
end

function free_ext_forces(N)
    """ Generates a trivial external force function """
    f_ext = zeros((N+1, 3))
    function f1(t, x, Q)
        return f_ext
    end
    return f1
end

function free_ext_torques()
    """ Generates a trivial external force function """
    function f1(t, x, Q)
        N1 = size(x, 1)
        return zeros((N1-1, 3))
    end
    return f1
end

function free_ext_torques(N)
    """ Generates a trivial external force function """
    t_ext = zeros((N, 3))
    function f1(t, x, Q)
        return t_ext
    end
    return f1
end


function interpolate_rodstate(x::Matrix{Float64}, Q::Array{Float64, 3}, Ns::Int)
    """ Utility function to interpolate x, Q to Nsx3 and (Ns-1)x3x3 arrays. """

    x_interp = zeros(Ns+1, 3)
    for i = 1:3
        interp_x = linear_interpolation(LinRange(0., 1.0, size(x,1)), x[:,i])
        x_interp[:,i] = interp_x.(LinRange(0., 1.0, Ns+1))
    end
    Q_interp = zeros(Ns, 3, 3)
    for i = 1:3
        for j=1:3
            interp_Q = linear_interpolation(LinRange(0., 1.0, size(x,1)-1), Q[:,i,j])
            Q_interp[:,i,j] = interp_Q.(LinRange(0., 1.0, Ns))
        end
    end
    # normalize Q
    for i = 1:Ns-1
        Q_interp[i,1,:] .= Q_interp[i,1,:] / norm(Q_interp[i,1,:])
        Q_interp[i,2,:] .= Q_interp[i,2,:] / norm(Q_interp[i,2,:])
        Q_interp[i,3,:] .= Q_interp[i,3,:] / norm(Q_interp[i,3,:])
    end

    return x_interp, Q_interp
end

function run_simulation(xpos0::Matrix{Float64}, Qs0::Array{Float64, 3},
                            dt::Float64, Ns::Int, T::Float64, n_saves::Int,
                            rod0::RodElements, rs0::Vector{Float64},
                            bcs::Any,
                            ext_forces::Any, ext_torques::Any, ell_rates::Vector{Float64}, plastic_rate::Float64, dxmax::Float64
                            )

    N = size(Qs0,1)



    # store data for saving - use the Ns parameter for a fixed saving size.
    # this is also the shape of the ext_forces and torques function output!
    xsaves = zeros(n_saves, Ns+1, 3)
    Qsaves = zeros(n_saves, Ns, 3, 3)
    tau = T / n_saves
    # BUGFIX: number of inner steps must be an integer loop bound.
    n_tsave = round(Int, T / dt / n_saves)

    rod = rod0
    rs = copy(rs0)
    # interpolate start configuration for saving

    x_i, Q_i = interpolate_rodstate(xpos0, Qs0, Ns)
    xsaves[1, :, :] = x_i[:,:]
    Qsaves[1, :,:,:] = Q_i[:,:,:]
    t = 0.0
    xpos = xpos0

    # normalize initial Q
    Qs = Qs0
    for i = 1:N
        Qs[i,1,:] .= Qs[i,1,:] / norm(Qs[i,1,:])
        Qs[i,2,:] .= Qs[i,2,:] / norm(Qs[i,2,:])
        Qs[i,3,:] .= Qs[i,3,:] / norm(Qs[i,3,:])
    end
    # initialize reference state
    l0s = get_l0s(xpos0, Qs)
    u0 = get_u0(xpos0, Qs, l0s)
    # create ell_rates_interp
    ell_rates_interp = linear_interpolation(LinRange(0,1.0, length(l0s)),ell_rates)
    @showprogress for nt = 1:(n_saves-1)
        for step = 1:n_tsave
            t += dt
            xpos, Qs, rod, u0, l0s, rs = timestepRK4plastic_remesh(t, xpos, Qs, l0s, u0, rs, dt,
                                        rod, bcs, ext_forces, ext_torques, ell_rates_interp, plastic_rate, dxmax)

        end
        # interpolate and save
        x_i, Q_i = interpolate_rodstate(xpos, Qs, Ns)
        xsaves[nt+1, :, :] = x_i[:,:]
        Qsaves[nt+1, :,:,:] = Q_i[:,:,:]
    end
    println("Simulation done")
    return xsaves, Qsaves
end


function timestepHeun(t::Float64, xpos::Matrix{Float64}, Qs::Array{Float64, 3}, dt::Float64, rod::RodElements, bcs::Function, ext_forces::Function, ext_torques::Function)

    # Heun method - predictor-corrector
    # stage 1: predictor
    forces1, torques1 = internal_forces_and_torques(xpos, Qs, rod)
    forces1 = forces1 + ext_forces(t, xpos, Qs)
    torques1 = torques1 + ext_torques(t, xpos, Qs)
    xnext1, Qnext1 = pos_next(xpos, Qs, forces1, torques1, dt, rod.gamma_p, rod.gamma_t, rod.N)
    xnext1, Qnext1 = bcs(t, xnext1, Qnext1)

    # stage 2: corrector
    forces2, torques2 = internal_forces_and_torques(xnext1, Qnext1, rod)
    forces2 = forces2 + ext_forces(t+dt, xnext1, Qnext1)
    torques2 = torques2 + ext_torques(t+dt, xnext1, Qnext1)
    # result
    xnext2, Qnext2 = pos_next(xpos, Qs, forces1+forces2, torques1+torques2, dt/2.0, rod.gamma_p, rod.gamma_t, rod.N)
    xnext2, Qnext2 = bcs(t+dt, xnext2, Qnext2)
    return xnext2, Qnext2

end



function remesh(xpos::Matrix{Float64}, Qs::Array{Float64, 3},
                l0s::AbstractArray{Float64}, u0::Matrix{Float64}, rod::RodElements, toobig::BitVector)

    N = size(Qs,1) # there are N segments, N+1 vertices

    N_new = N + sum(toobig)

    xpos_new = zeros(N_new+1,3)
    Qs_new = zeros(N_new,3,3)
    l0s_new = zeros(N_new)
    u0_new = zeros(N_new-1, 3) # only on interior points!
    # assumption: Stretching and shearing moduli are large and uniform.
    # so we only need to update the bend and twist moduli.
    B_new = zeros(N_new)
    Bt_new = zeros(N_new)
    #rs_new = zeros(N_new)
    idx_new = 0
    u0_m = zeros(3)
    # initial point doesn't move
    for idx = 1:N
        idx_new +=1
        if toobig[idx]
            #update current - this one has same position as current
            xpos_new[idx_new,:] = xpos[idx,:]
            Qs_new[idx_new,:,:] = Qs[idx,:,:]
            #rs_new[idx_new] = rs[idx]
            if idx > 1
                u0_new[idx_new-1,:] = u0[idx-1,:]  # voro adjustment factor?
            end
            l0s_new[idx_new] = 0.5*l0s[idx] # split in half
            # BUGFIX: the original (first half of the split) segment must keep the
            # parent's bend/twist moduli, otherwise it is left at zero.
            B_new[idx_new] = rod.B[idx]
            Bt_new[idx_new] = rod.Bt[idx]

            idx_new +=1
            # make the new segment: it is midway through both, colinear with Q,
            # and initially straight with same radius as the parent edge
            xpos_new[idx_new,:] = 0.5*(xpos[idx,:]+xpos[idx+1,:])
            Qs_new[idx_new,:,:] = Qs[idx,:,:]
            #rs_new[idx_new] = rs[idx]
            u0_new[idx_new-1,:] = [0., 0., 0.]
            l0s_new[idx_new] = 0.5*l0s[idx]
            #B_new[idx_new] = 0.5*(rod.B[idx]+rod.B[idx+1])
            #Bt_new[idx_new] = 0.5*(rod.Bt[idx]+rod.Bt[idx+1])
            B_new[idx_new] = rod.B[idx]
            Bt_new[idx_new] = rod.Bt[idx]

        else
            xpos_new[idx_new,:] = xpos[idx,:]
            Qs_new[idx_new,:,:] = Qs[idx,:,:]
            #rs_new[idx_new] = rs[idx]
            l0s_new[idx_new] = l0s[idx]
            # BUGFIX: carry over moduli for unsplit segments too.
            B_new[idx_new] = rod.B[idx]
            Bt_new[idx_new] = rod.Bt[idx]
            if idx > 1
                u0_new[idx_new-1,:] = u0[idx-1,:]  # voro adjustment factor?
            end
        end
    end
    # take care of the end point
    xpos_new[end,:] = xpos[end,:] # this one doesn't move.

    return xpos_new, Qs_new, l0s_new, u0_new, B_new, Bt_new
end

function split_weights(w::Vector{Float64}, toobig::BitVector)
    """Fraction of an initial edge carried by each edge after remesh splits the edges
    flagged in toobig (same ordering as remesh: each split edge -> two halves in place)"""
    w_new = Float64[]
    for k in eachindex(w)
        if toobig[k]
            push!(w_new, 0.5 * w[k], 0.5 * w[k])
        else
            push!(w_new, w[k])
        end
    end
    return w_new
end


function timestepRK4plastic_remesh(t::Float64, xpos::Matrix{Float64}, Qs::Array{Float64, 3},
                            l0s::AbstractArray{Float64}, u0::Matrix{Float64}, dt::Float64,
                            rod::RodElements, bcs::Function, neighbour_walls::Vector{Vector{Plane}},
                            alpha::Float64, neighbour_list::Vector{Vector{Int}}, cradii::Vector{Float64},
                            ext_forces::Function, ext_torques::Function,
                            ell_rates, plastic_rate::Float64, dxmax::Float64)

    ## Runge-Kutta order 4
    # stage 1
    forces1, torques1 = internal_forces_and_torques(xpos, Qs, rod, u0, l0s)
    forces1 += ext_forces(xpos, neighbour_walls)
    forces1 += compute_contact_force_clamped(xpos, neighbour_list, cradii, alpha)
    torques1 += ext_torques(t, xpos, Qs)
    xnext1, Qnext1 = pos_next(xpos, Qs, forces1, torques1, dt/2.0, rod.gamma_p, rod.gamma_t, rod.N)
    xnext1, Qnext1 = bcs(t, xnext1, Qnext1)
    # stage 2
    forces2, torques2 = internal_forces_and_torques(xnext1, Qnext1, rod, u0, l0s) # should we use updated u0, ls0?
    forces2 += ext_forces(xnext1, neighbour_walls)
    forces2 += compute_contact_force_clamped(xnext1, neighbour_list, cradii, alpha)
    torques2 += ext_torques(t+dt/2.0, xnext1, Qnext1)
    xnext2, Qnext2 = pos_next(xpos, Qs, forces2, torques2, dt/2.0, rod.gamma_p, rod.gamma_t, rod.N)
    xnext2, Qnext2 = bcs(t+dt/2.0, xnext2, Qnext2)
    # stage 3
    forces3, torques3 = internal_forces_and_torques(xnext2, Qnext2, rod, u0, l0s)
    forces3 += ext_forces(xnext2, neighbour_walls)
    forces3 += compute_contact_force_clamped(xnext2, neighbour_list, cradii, alpha)
    torques3 += ext_torques(t+dt/2.0, xnext2, Qnext2)
    xnext3, Qnext3 = pos_next(xpos, Qs, forces3, torques3, dt, rod.gamma_p, rod.gamma_t, rod.N)
    xnext3, Qnext3 = bcs(t+dt/2.0, xnext3, Qnext3)
    # stage 4
    forces4, torques4 = internal_forces_and_torques(xnext3, Qnext3, rod, u0, l0s)
    forces4 += ext_forces(xnext3, neighbour_walls)
    forces4 += compute_contact_force_clamped(xnext3, neighbour_list, cradii, alpha)
    torques4 += ext_torques(t+dt, xnext3, Qnext3)
    # result
    xnext_RK4, Qnext_RK4 = pos_next(xpos, Qs, forces1 + 2.0*forces2 + 2.0*forces3 + forces4,
                                            torques1 + 2.0*torques2 + 2.0*torques3 + torques4, dt/6.0, rod.gamma_p, rod.gamma_t, rod.N)
    xnext_RK4, Qnext_RK4 = bcs(t+dt, xnext_RK4, Qnext_RK4)
    ## plasticity and growth
    #update u0, l0s:
    #ell_rates = ell_rates_interp.(LinRange(0,1,length(l0s)))
    #l0s_new = (1.0 .+ ell_rates * dt ) .* l0s
    l0s_new = l0s .+ ell_rates * dt
    #rs_new = rs ./ (1.0 .+ ell_rates * dt)
    u0_new = get_u0_plastic(xnext_RK4, Qnext_RK4, l0s_new,  u0, dt, plastic_rate)

    return xnext_RK4, Qnext_RK4, u0_new, l0s_new

end



"""
    fix_total_growth(ell_rates, w_seg, total_growth, t)

Optional normalization of the growth (keyword `total_growth` of run_simulation_plastic_remeshed):
`total_growth(t)` is the prescribed rate of change of the total rest length. The edge elongation
rates are scaled so that they sum to it, keeping their profile along the rod. If the summed rate is
not clearly positive (< half the target), a uniform rate per initial edge (∝ w_seg) is added instead.
`total_growth = nothing` (default) leaves the rates unchanged.
"""
fix_total_growth(ell_rates, w_seg, ::Nothing, t) = ell_rates
function fix_total_growth(ell_rates, w_seg, total_growth, t)
    target = total_growth(t)
    tot = sum(ell_rates)
    if tot > 0.5 * target && tot > 0
        return ell_rates .* (target / tot)
    else
        return ell_rates .+ w_seg .* ((target - tot) / sum(w_seg))
    end
end

function run_simulation_plastic_remeshed(xpos0::Matrix{Float64}, Qs0::Array{Float64, 3},
    dt::Float64, Ns::Int, T::Float64, n_saves::Int,
    rod0::RodElements, S::Float64, St::Float64,
    bcs::Any, circum_radii, alpha::Float64,
    walls::Vector{Plane}, cutoff::Float64, ext_forces::Any, ext_torques::Any, growth_rates, plastic_rate::Float64, dxmax::Float64; elastic=false, proba=0.0,
    total_growth=nothing
    )

    N = size(Qs0,1)



    # store data for saving - use the Ns parameter for a fixed saving size.
    # this is also the shape of the ext_forces and torques function output!
    xsaves = zeros(n_saves, Ns+1, 3)
    Qsaves = zeros(n_saves, Ns, 3, 3)
    tau = T / n_saves
    # BUGFIX: number of inner steps must be an integer loop bound.
    n_tsave = round(Int, T / dt / n_saves)

    rod = rod0
    # interpolate start configuration for saving

    x_i, Q_i = interpolate_rodstate(xpos0, Qs0, Ns)
    xsaves[1, :, :] = x_i[:,:]
    Qsaves[1, :,:,:] = Q_i[:,:,:]
    t = 0.0
    xpos = xpos0

    # normalize initial Q
    Qs = Qs0
    for i = 1:N
        Qs[i,1,:] .= Qs[i,1,:] / norm(Qs[i,1,:])
        Qs[i,2,:] .= Qs[i,2,:] / norm(Qs[i,2,:])
        Qs[i,3,:] .= Qs[i,3,:] / norm(Qs[i,3,:])
    end
    # initialize reference state

    l0s = get_l0s(xpos0, Qs)
    if proba == 0.
        if elastic
            u0 = zeros(N-1,3)
        else
            u0 = get_u0(xpos0, Qs, l0s)
        end
    else
        u0 = zeros(N-1,3)
        u_r = get_u0(xpos0, Qs, l0s)
        for k = 1:(N-1)
            u0[k, :] .= u0[k, :] + proba[k] * u_r[k,:]
        end
    end

    s_coords = get_fractional_arclength(xpos)
    # initialize walls list
    neighbour_walls = make_neighbour_walls(xpos, walls, cutoff)

    # BUGFIX (2026-10-19): growth_rates(t, s) is the elongation rate of an edge of the
    # initial discretization. Each edge carries the fraction w_seg of an initial edge
    # (halved at every split by remesh), and grows at w_seg * growth_rates. Previously both
    # children of a split edge grew at the parent's full rate, doubling the local growth.
    w_seg = ones(N)


    # 1st time point s_coords = get_fractional_arclength(xpos)
            itp_points = [(t, s) for s = s_coords]
            cradii = [circum_radii(points...) for points = itp_points]

            neighbour_list = make_neighbour_points(xpos, cradii)
            ell_rates = fix_total_growth(w_seg .* [growth_rates(points...) for points = itp_points][2:end], w_seg, total_growth, t)
            xpos, Qs, u0, l0s = timestepRK4plastic_remesh(t, xpos, Qs, l0s, u0, dt, rod, bcs, neighbour_walls,
                                                            alpha, neighbour_list, cradii,
                                                             ext_forces, ext_torques, ell_rates, plastic_rate, dxmax)
    @showprogress for nt = 1:(n_saves-1)
        for step = 1:n_tsave
            t += dt

            #ell_rates_interp = linear_interpolation(s_coords, vec(ell_rates(t)); extrapolation_bc=Flat())
            #ell_rates_interp = linear_interpolation(LinRange(0,1.0, Ns), vec(ell_rates(t)))
            if step % 5 == 0
                # remesh and update walls
                #println("Checking for remeshing at time ", t)
                toobig = l0s .> dxmax
                if any(toobig)
                    #println("Remeshing at time ", t, " with ", sum(toobig), " segments to be split.")
                    # insert new elements for each segment that got too big.
                    xpos, Qs, l0s, u0, Bvec, Btvec  = remesh(xpos, Qs, l0s, u0, rod, toobig)
                    w_seg = split_weights(w_seg, toobig)
                    # assemble new rod with new values of radius, and right number of segments after remeshing
                    rod = make_Rod(xpos, S, St, Bvec, Btvec, rod.gamma_p, rod.gamma_t)
                end
                neighbour_walls = make_neighbour_walls(xpos, walls, cutoff)
                s_coords = get_fractional_arclength(xpos)
                itp_points = [(t, s) for s = s_coords]
                ell_rates = fix_total_growth(w_seg .* [growth_rates(points...) for points = itp_points][2:end], w_seg, total_growth, t)
                cradii = [circum_radii(points...) for points = itp_points]

                neighbour_list = make_neighbour_points(xpos, cradii)
            else
                s_coords = get_fractional_arclength(xpos)
                itp_points = [(t, s) for s = s_coords]
                ell_rates = fix_total_growth(w_seg .* [growth_rates(points...) for points = itp_points][2:end], w_seg, total_growth, t)
                cradii = [circum_radii(points...) for points = itp_points]
            end
            xpos, Qs, u0, l0s = timestepRK4plastic_remesh(t, xpos, Qs, l0s, u0, dt, rod, bcs, neighbour_walls,
            alpha, neighbour_list, cradii, ext_forces, ext_torques, ell_rates, plastic_rate, dxmax)

        end
            # interpolate and save
            x_i, Q_i = interpolate_rodstate(xpos, Qs, Ns)
            xsaves[nt+1, :, :] = x_i[:,:]
            Qsaves[nt+1, :,:,:] = Q_i[:,:,:]
        end
    println("Simulation done")
    return xsaves, Qsaves
end
