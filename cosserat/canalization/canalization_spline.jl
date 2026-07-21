## Slow timescale simulations
#

using Distributed

@everywhere begin import Pkg; Pkg.activate("../src/rods")
    Pkg.instantiate();
end
@everywhere begin # imports
    include("../src/rod_contact.jl")
    include("../src/walls.jl")
    include("../src/contact.jl")
    using ProgressMeter
    using LinearAlgebra
    using Random
    using HDF5   
    using MAT 
    using Statistics
    using Interpolations
end


@everywhere begin # parameter definitions
    
    path = "./SRVbox2f_16_notail/"
    key = ".mat"
    mat_files = filter(x -> occursin(key, x), readdir(path))

    num_files = length(mat_files)
    

end

println("---- Starting ----")
status = @showprogress pmap(1:num_files) do n
    
    # here run simulations with the corresponding parameter values.

    # pick the initial condition with the n

    #load the data
    vars = matread(path * mat_files[n])
    tau = vars["tau"]
    data = vars["x0"]
    x_c = 1
    t0 = 1
    x = data[x_c:2:end-1,:] # the first timepoints are worse...

    Q = vars["Q_i"]
    println("Size Q:", size(Q))
    Qstart = Q[t0, x_c:2:end-1, :, :]
    # adapt Qstart to x
    tan0 = x[2:end,:]-x[1:end-1,:]
    tan0 = tan0 ./ sqrt.(sum(abs2, tan0, dims=2))
    Qstart[:,1,:] .= Qstart[:,1,:] - dot.(eachrow(Qstart[:,1,:]), eachrow(tan0)) .* tan0
    Qstart[:,2,:] .= Qstart[:,2,:] - dot.(eachrow(Qstart[:,2,:]), eachrow(tan0))  .* tan0
    Qstart[:,3,:] .= tan0

    Qstart[:,1,:] = Qstart[:,1,:] ./ sqrt.(sum(abs2, Qstart[:,1,:], dims=2))
    Qstart[:,2,:] = Qstart[:,2,:] ./ sqrt.(sum(abs2, Qstart[:,2,:], dims=2))
    Qstart[:,3,:] = Qstart[:,3,:] ./ sqrt.(sum(abs2, Qstart[:,3,:], dims=2))

    
    num_t = size(Q,1)
    data = vars["x_i"]
    x_itp = LinearInterpolation(LinRange(0, 1, num_t-t0+1), [data[t,x_c:2:end-1,:] for t = t0:num_t])


    Q1_itp = LinearInterpolation(LinRange(0, 1, num_t-t0+1), [Q[t,1,:,:] for t = t0:num_t])
    Qend_itp = LinearInterpolation(LinRange(0, 1, num_t-t0+1), [Q[t,end,:,:] for t = t0:num_t])
    
    # boundary polyhedron definition
    normals = vars["n"]
    ds = vars["d"]
    
    # compute growth rate
    num_s = size(x, 1)
    segment_lengths = [sqrt.(sum(abs2,data[t,x_c:2:end-1,:][2:end,:] -data[t,x_c:2:end-1,:][1:end-1,:];dims=2)) for t = t0:num_t]
    arclengths = [cumsum([0.0; segment_lengths[t]], dims=1) for t = 1:length(segment_lengths)]
    length_t = [arclengths[t][end] for t = 1:length(arclengths)]
    fractional = [arclengths[t] ./ length_t[t] for t = 1:length(arclengths)]
    # smooth the growth rate - this is also reducing the magnitude.
    for k = 1:4
        segment_lengths[2:end-1] += 0.5*(segment_lengths[1:end-2]+segment_lengths[3:end] - 2* segment_lengths[2:end-1])
       end
    delta_lengths = [segment_lengths[t+1] - segment_lengths[t] for t =1:length(segment_lengths)-1]
    delta_lengths_scoord = [LinearInterpolation(vec(fractional[t][1:end-1]), vec(delta_lengths[t]); extrapolation_bc=Flat()) for t =1:length(delta_lengths)]

    growth_arr = zeros(num_t-t0, num_s)
    s_i = [s for s = LinRange(0, 1, num_s)]
    for t = 1:(num_t-t0)
        growth_arr[t, :] = delta_lengths_scoord[t](s_i)
    end

    dt_growth = 1.0/(num_t-t0)

    ## interpolate circumradii
    c_radii = vars["circum_radii"][t0:end, x_c:2:end-1]
    cradii_arr = zeros(num_t-t0, num_s)
    c_radii_scoord = [LinearInterpolation(vec(fractional[t][1:end]), vec(c_radii[t, 1:end]); extrapolation_bc=Flat()) for t =1:length(delta_lengths)]
    for t = 1:(num_t-t0)
        cradii_arr[t, :] = 1.0*c_radii_scoord[t](s_i)
    end
    
    T = 1.0
    twist_ratio = 2.0
    plastic_rate = 10.0

    #twist_ratio = 1.0

    asym = 1.0
    x[:,2] = asym * x[:,2]
    max_x = maximum(x[:,1])
    
    
   N = size(x,1)-1
   Ns = N

   L = sum([norm(x[i+1,:]-x[i,:]) for i=1:size(x,1)-1])
    
   

   x0 = copy(x)
   x0 = x0 / L
    
    dxmax = 2.0*maximum([norm(x0[i+1,:]-x0[i,:]) for i=1:size(x,1)-1])
    
    walls = make_walls(normals, ds/L) # need to rescale d by length unit
    # apply perturbation
    
    
    # confinement potential specification
    
    r0_um = 0.5*80/L

   nu = 0.3
   S = 1e3 #10.0 #1e3
   St = 1e3
   B = 1.0 #[B_ave - B_diff*tanh(20*(i/N-0.75)) for i = 1:N] # linear tapering of bending stiffness
   Bt = twist_ratio * B
   alpha = 0*3e2 # self-contact parameter

   gamma_p = 1.0
   gamma_t = 1.0
    
   dx = minimum([norm(x0[i+1,:]-x0[i,:]) for i=1:size(x0,1)-1])
   rod = make_Rod(x0, S, St, B, Bt, gamma_p, gamma_t) # adapt 
   x_start = copy(x0)
   #Qstart = get_adaptedQ(x0, 0.)
   q_0 = copy(Qstart[1,:,:])
   q_1 = copy(Qstart[end,:,:])

    #T = values_vec[n] # set above
    #dt = 1e-3* dx #
    dt = 3e-4* dx
    n_saves = 200
    
    ## make growth rate interpolant
    nodes = ([t for t = LinRange(0, T, num_t-t0)], [s for s = LinRange(0, 1, num_s)])
    growth_rate = extrapolate(interpolate(nodes, growth_arr / (L*dt_growth*T), Gridded(Linear())), Flat())
    circum_radii = extrapolate(interpolate(nodes, cradii_arr / L, Gridded(Linear())), Flat())

    
    println("dt = $(dt)")
    
   function bcs(t, x, Q) # no-twist, end clamped bcs
        xp = copy(x) 
        Qp = copy(Q)
        Qp[1,:,:] = copy(Q1_itp(t/T)) #RodriguesMat(Qstart[1,:,:], omega0*t, Qstart[1,:,3])
        Qp[end,:,:] = copy(Qend_itp(t/T)) #Rodrigues
        xp[1,:] = x_itp(t/T)[1,:]/L #copy(x_start[1,:]) # initial anterior point
        #xp[1,1] += (t/T) * j_x # jogging of x-axis of anterior point
        #xp[1,2] += (t/T) * j_y # jogging of y-axis of anterior point
        #xp[1,3] += (t/T) * j_z # jogging of z-axis of anterior point
        xp[end,:] = x_itp(t/T)[end,:]/L #copy(x_start[end,:])
        return xp, Qp
   end

    

    #tau = zeros(size(x,1)-1,3)
    xi = collect(LinRange(-1, 1, N))
    
    # save parameters
    params = [T, plastic_rate, S, St, alpha, twist_ratio]

    function ext_torques(t::Float64,x::Matrix{Float64}, Q::Array{Float64, 3})
        return zeros(size(Q,1), 3)
    end

    function ext_forces(x, neighbour_walls) # confinement forces
        f_ext = zeros(size(x,1),3)
        N = size(x,1)-1
        for i = 2:N
            force_orient, dist = get_force_components(x[i,:], neighbour_walls[i])
            f_ext[i,:] = 10000 * dist * (1 + 00*dist) * force_orient
        end
        return f_ext
    end
    
    cutoff = 0.1
    
    xsaves, Qsaves = run_simulation_plastic_remeshed(x_start, Qstart, dt, Ns, T, n_saves, 
        rod, S, St, bcs, circum_radii, alpha, walls, cutoff, ext_forces, ext_torques, growth_rate, plastic_rate, dxmax)
    
    #out_filename = "./plas2/cutoff01plastic10_T$(T).mat"
    folder_name = "./20260402_canalization_box2f_nocontact_notail_s6s1"
    mkpath(folder_name)
    out_filename = folder_name * "/cutoff01contact_n$(n)_tau$(tau).mat"
    matwrite(out_filename, Dict("xsaves" => xsaves, "Qsaves" => Qsaves, "params" => params))
end
println("~success~")