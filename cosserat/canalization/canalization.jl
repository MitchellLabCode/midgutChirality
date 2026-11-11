## Canalization sweep: rods grown from initial shapes interpolated between two imaged midguts
#
# One run per input file <CANAL_INPUT>/itp_<k>_tau<tau>.mat (tau = 0: sample 6, WT; tau = 1: sample 3,
# byn>Myo1C; alpha = 2 tau - 1). Each file holds the interpolated time series x_i / Q_i, which sets the
# clamped ends and the growth-rate field (profile along the rod and timing), the start shape x0, the
# circumradii (self-contact) and the confining walls n, d (SRVbox3f_s6s3: one polyhedron containing every
# interpolated shape at every frame).
#
# By default the total rest-length growth rate is prescribed (keyword total_growth), so that every run
# reaches the same final length L/L0: the growth of the input series dips in the middle of the sweep
# (an effect of the SRV interpolation), which would otherwise tie the final length to alpha.
#
#   CANAL_GROWTH = common (default)  every run grows to R = mean of the input growth of the two
#                                    measured ends (the files nearest tau = 0 and tau = 1)
#                = <number>          every run grows to R = that number
#                = none              no normalization (control)
#   CANAL_INPUT  = input folder (default SRVbox3f_s6s3: sample 6 <-> sample 3)
#   CANAL_START  = xi (default: start from the first frame of the input series x_i)
#                | x0 (start from the x0 field of the input file, as the April 2026 runs; the clamped ends
#                  and growth still follow x_i, lengths are in units of the length of x0)
#   CANAL_ALPHA  = self-contact stiffness (default 300); other values add _a<alpha> to the folder name
#   CANAL_DATE   = date stamp of the output folder (default 20261110)
#   CANAL_SOLVER = fast (default: ../prod/fastprod.jl, ~10x faster, same results to round-off)
#                | ref  (run_simulation_plastic_remeshed of ../src/rod_contact2.jl)
#                | both (runs both, prints max |dx|, |dQ| and the timings, saves the fast result)
#   CANAL_ONLY   = comma-separated run indices n (default: all files)
#   CANAL_T, CANAL_NSAVES = simulation time and number of saves (default 3.0, 200)
# Validation, e.g.: CANAL_SOLVER=both CANAL_T=0.1 CANAL_NSAVES=10 CANAL_ONLY=1,25,50 CANAL_DATE=validate \
#                   julia -p 3 canalization.jl   (or: sbatch validate_fast.sh)
#
# The total rate keeps the time course of the run's own input series: dL/dt(t) ∝ the summed length
# change of its frames, scaled so that ∫ dL/dt dt = R - 1 (L0 = 1).
#

using Distributed

@everywhere begin import Pkg; Pkg.activate("../src/rods")
    Pkg.instantiate();
end
@everywhere begin # imports
    include("../src/rod_contact2.jl")
    include("../src/walls.jl")
    include("../src/contact.jl")
    include("../src/fastkernel.jl")              # fast solver: StaticArrays kernel + run_fastprod
    include("../prod/fastprod.jl")
    using ProgressMeter
    using LinearAlgebra
    using Random
    using HDF5   
    using MAT 
    using Statistics
    using Interpolations
end


@everywhere begin # parameter definitions
    
    # input series: SRVbox3f_s6s3 = sample 6 (tau = 0) <-> sample 3 (tau = 1)
    canal_input = get(ENV, "CANAL_INPUT", "SRVbox3f_s6s3")
    box_tag  = replace(split(canal_input, "_")[1], "SRV" => "")           # box3f
    pair_suffix = occursin("_", canal_input) ? String(split(canal_input, "_")[2]) : ""
    pair_tag = startswith(pair_suffix, "s") ? pair_suffix :          # e.g. SRVbox3f_s6s3
               get(Dict("16" => "s6s1", "56" => "s6s5", "" => "s6s3"), pair_suffix, "s" * pair_suffix)
    path = "./$(canal_input)/"
    key = ".mat"
    mat_files = filter(x -> occursin(key, x), readdir(path))

    num_files = length(mat_files)

    growth_mode = get(ENV, "CANAL_GROWTH", "common")
    start_mode  = get(ENV, "CANAL_START", "xi")      # xi: first frame of x_i;  x0: the x0 field (as the April runs)
    contact_alpha = parse(Float64, get(ENV, "CANAL_ALPHA", "300"))   # self-contact stiffness (canalization.jl: 3e2)
    canal_date  = get(ENV, "CANAL_DATE", "20261110")
    solver      = get(ENV, "CANAL_SOLVER", "fast")
    solver in ("fast", "ref", "both") || error("CANAL_SOLVER = fast | ref | both")
    T_sim       = parse(Float64, get(ENV, "CANAL_T", "3.0"))
    n_saves_sim = parse(Int, get(ENV, "CANAL_NSAVES", "200"))
    run_list    = haskey(ENV, "CANAL_ONLY") ? parse.(Int, split(ENV["CANAL_ONLY"], ",")) : collect(1:num_files)
    # growth L_end / L_0 of an input series (centerline subsampled as below)
    function input_growth(file)
        xi = matread(path * file)["x_i"][:, 1:2:end-1, :]
        lens = [sum(norm.(eachrow(xi[t, 2:end, :] - xi[t, 1:end-1, :]))) for t = 1:size(xi, 1)]
        return lens[end] / lens[1]
    end
    if growth_mode == "common"
        taus = [matread(path * f)["tau"] for f in mat_files]
        i0, i1 = argmin(abs.(taus .- 0.0)), argmin(abs.(taus .- 1.0))
        R_target = 0.5 * (input_growth(mat_files[i0]) + input_growth(mat_files[i1]))
    elseif growth_mode == "none"
        R_target = NaN
    else
        R_target = parse(Float64, growth_mode)
    end
    folder_tag = (growth_mode == "none" ? "_ctrl" : (growth_mode == "common" ? "_Lnorm" : "_L$(growth_mode)")) *
                 (start_mode == "x0" ? "_x0" : "") *
                 (contact_alpha == 300 ? "" : "_a$(contact_alpha == round(contact_alpha) ? Int(contact_alpha) : contact_alpha)")
    folder_name = "./$(canal_date)_canalization_$(box_tag)_$(pair_tag)$(folder_tag)"


end

println("---- Starting ----")
status = @showprogress pmap(run_list) do n
    
    # here run simulations with the corresponding parameter values.

    # pick the initial condition with the n

    #load the data
    vars = matread(path * mat_files[n])
    tau = vars["tau"]
    data = vars["x_i"]
    x_c = 1
    t0 = 1
    x = data[t0, x_c:2:end-1,:] # the first timepoints are worse...
    if start_mode == "x0"       # start from the x0 field instead (the length unit L below is then its length)
        x = Matrix{Float64}(vars["x0"][x_c:2:end-1, :])
    end

    Q = vars["Q_i"]
    Qstart = Q[t0, x_c:2:end-1, :, :]
    
    num_t = size(data,1)
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
    
    T = T_sim
    twist_ratio = 1.0
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
   alpha = contact_alpha # self-contact stiffness (default 3e2)

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
    n_saves = n_saves_sim
    
    ## make growth rate interpolant
    nodes = ([t for t = LinRange(0, T, num_t-t0)], [s for s = LinRange(0, 1, num_s)])
    growth_rate = extrapolate(interpolate(nodes, growth_arr / (L*dt_growth*T), Gridded(Linear())), Flat())
    circum_radii = extrapolate(interpolate(nodes, cradii_arr / L, Gridded(Linear())), Flat())

    ## prescribed total growth rate: the summed length change of the input frames (its time
    ## course), scaled so that the rest length grows from 1 to R_target
    if growth_mode == "none"
        total_growth = nothing
        R_input = NaN
    else
        D = [sum(delta_lengths[k]) for k = 1:(num_t-t0)] / (L*dt_growth*T)   # as growth_arr
        D_itp = extrapolate(interpolate((nodes[1],), D, Gridded(Linear())), Flat())
        tt = LinRange(0, T, 20001)
        I = sum(D_itp.(tt)) * (tt[2] - tt[1]) - 0.5 * (D_itp(0.0) + D_itp(T)) * (tt[2] - tt[1])
        R_input = 1 + I
        total_growth = t -> D_itp(t) * (R_target - 1) / I
    end
    println("run $(n): tau = $(tau), input growth $(round(R_input, digits=3)) -> target $(R_target)")

    
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
    
    function bc!(t, xv, Qv) # the same clamped ends, in place on Vector{V3} / Vector{M3} (fastprod)
        Qv[1] = M3(Q1_itp(t/T))
        Qv[end] = M3(Qend_itp(t/T))
        xi_t = x_itp(t/T)
        xv[1] = V3(xi_t[1,:]/L)
        xv[end] = V3(xi_t[end,:]/L)
        return nothing
    end
    run_ref() = run_simulation_plastic_remeshed(copy(x_start), copy(Qstart), dt, Ns, T, n_saves,
        rod, S, St, bcs, circum_radii, alpha, walls, cutoff, ext_forces, ext_torques, growth_rate, plastic_rate, dxmax;
        total_growth = total_growth)
    run_fast() = run_fastprod(copy(x_start), copy(Qstart), dt, Ns, T, n_saves, S, St, B, Bt, gamma_p, gamma_t,
        bc!, circum_radii, alpha, walls, cutoff, 10000.0, growth_rate, plastic_rate, dxmax;
        material_u0 = true, total_growth = total_growth)     # walls: 10000 * depth, as ext_forces
    t_start = time()
    if solver == "ref"
        xsaves, Qsaves = run_ref()
    else
        xsaves, Qsaves = run_fast()
    end
    elapsed = time() - t_start
    if solver == "both"
        t_ref = time(); xr, Qr = run_ref(); t_ref = time() - t_ref
        println("VALIDATE run $(n) (tau = $(round(tau, digits=4))): max |dx| = $(maximum(abs.(xr - xsaves))), " *
                "max |dQ| = $(maximum(abs.(Qr - Qsaves))) over $(n_saves) saves; " *
                "time ref $(round(t_ref, digits=1)) s, fast $(round(elapsed, digits=1)) s ($(round(t_ref / elapsed, digits=1))x)")
    end
    
    #out_filename = "./plas2/cutoff01plastic10_T$(T).mat"
    mkpath(folder_name)
    out_filename = folder_name * "/cutoff01contact_n$(n)_tau$(tau).mat"
    matwrite(out_filename, Dict("xsaves" => xsaves, "Qsaves" => Qsaves, "params" => params,
                                "growth_mode" => growth_mode, "input" => canal_input, "start" => start_mode, "growth_target" => R_target, "growth_input" => R_input,
                                "solver" => solver == "ref" ? "ref" : "fast", "elapsed_s" => elapsed))
end
println("~success~")