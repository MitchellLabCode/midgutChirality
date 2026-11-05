## Shared setup for the production runs, with the V6 reference data (refIV6_sample*_walls3_cball.mat)
# as inputs: sample registry, parameter grid, per-run parameters (prepare_run) and the clamped,
# data-driven ends. The initial centerlines x0_s*V5_t*.mat (in this folder) are the V5 ones.
# Included by QTpTw_batch_samples_V6.jl and validate_fastprod.jl.

x_c = 1

# ---- sample registry --------------------------------------------------------
# ref  : reference .mat file (holds Q, x, circum_radii, n, d)
# x0   : initial centerline file (holds "x0")
# t0   : start time index into the reference data
# tag  : short label used in the output folder name
sample_configs = Dict(
    1 => (ref = "../data/refIV6_sample1_walls3_cball.mat", x0 = "./x0_s1V5_t38_ant0.mat",     t0 = 38, tag = "s1V6"),
    2 => (ref = "../data/refIV6_sample2_walls3_cball.mat", x0 = "./x0_s2V5_t25_normal.mat",    t0 = 25, tag = "s2V6"),
    3 => (ref = "../data/refIV6_sample3_walls3_cball.mat", x0 = "./x0_s3V5_t45_f3.mat",        t0 = 45, tag = "s3V6"),
    4 => (ref = "../data/refIV6_sample4_walls3_cball.mat", x0 = "./x0_s4V5_t36_normal.mat",    t0 = 36, tag = "s4V6"),
    5 => (ref = "../data/refIV6_sample5_walls3_cball.mat", x0 = "./x0_s5V5_t25_normal.mat",    t0 = 25, tag = "s5V6"),
    6 => (ref = "../data/refIV6_sample6_walls3_cball.mat", x0 = "./x0_s6V5_t54_ant.mat",       t0 = 54, tag = "s6V6"),
    7 => (ref = "../data/refIV6_sample7_walls3_cball.mat", x0 = "./x0_s7V5_t15_normal.mat", t0 = 15, tag = "s7V6"),
    71 => (ref = "../data/refIV6_sample7_walls3_cball.mat", x0 = "./x0_s7V5_t20_normal.mat", t0 = 20, tag = "s7V6"),
    8 => (ref = "../data/refIV6_sample8_walls3_cball.mat", x0 = "./x0_s8V5_t37_normal.mat", t0 = 37, tag = "s8V6"),
    9 => (ref = "../data/refIV6_sample9_walls3_cball.mat", x0 = "./x0_s9V5_t45_normal.mat", t0 = 45, tag = "s9V6"),
    10 => (ref = "../data/refIV6_sample10_walls3_cball.mat", x0 = "./x0_s10V5_t45_normal.mat", t0 = 45, tag = "s10V6"),
    # V6 addition: sample 10 from t0 = 35, the start used for prod's selected/ runs
    101 => (ref = "../data/refIV6_sample10_walls3_cball.mat", x0 = "./x0_s10V5_t35_normal.mat", t0 = 35, tag = "s10V6")
)

# samples used when no ids are given on the command line
default_samples = [1, 2, 3, 4, 5, 6, 7, 71, 8, 9, 10, 101]

# date stamp used in the output folder name
date_str = "20261019"

# ---- parameter sweep (shared across samples) --------------------------------
# Fine log-spaced grid.
#Ts                = [1.0, 3.0] # 10 .^ collect(-1:0.2:log10(30))
#twist_ratios      = [0.; 10 .^ collect(-1:0.1:1)]
#plastic_rates_arr = [0.; 10 .^ collect(-2:0.3:3)]

# Coarse alternative (matching ./Qtwp_contact_batch.jl) -- swap in if the fine
# grid above is too expensive with self-contact enabled:
Ts                = [1.0, 3.0]
twist_ratios      = [0.,  0.3, 1.0, 3.0]
plastic_rates_arr = [0.,  0.1, 1.0, 10.0, 100.0, 1000.0]

values_vec = [(T, tw, pr) for T in Ts for tw in twist_ratios for pr in plastic_rates_arr]

# ---- per-sample setup -------------------------------------------------------
# Builds everything the simulation closure needs for a given sample and returns
# it as a NamedTuple. Mirrors the parameter block of Qtwp_contact_batch.jl.
function setup_sample(sample_id::Int)
    cfg = sample_configs[sample_id]

    vars = matread(cfg.ref)
    data = vars["x"]
    t0 = cfg.t0

    data_0 = matread(cfg.x0)["x0"]
    x = data_0[x_c:2:end, :]

    Q = vars["Q"]
    Qstart = Q[t0, x_c:2:end-1, :, :]
    # adapt Qstart to x
    tan0 = x[2:end, :] - x[1:end-1, :]
    tan0 = tan0 ./ sqrt.(sum(abs2, tan0, dims=2))
    Qstart[:, 1, :] .= Qstart[:, 1, :] - dot.(eachrow(Qstart[:, 1, :]), eachrow(tan0)) .* tan0
    Qstart[:, 2, :] .= Qstart[:, 2, :] - dot.(eachrow(Qstart[:, 2, :]), eachrow(tan0)) .* tan0
    Qstart[:, 3, :] .= tan0

    Qstart[:, 1, :] = Qstart[:, 1, :] ./ sqrt.(sum(abs2, Qstart[:, 1, :], dims=2))
    Qstart[:, 2, :] = Qstart[:, 2, :] ./ sqrt.(sum(abs2, Qstart[:, 2, :], dims=2))
    Qstart[:, 3, :] = Qstart[:, 3, :] ./ sqrt.(sum(abs2, Qstart[:, 3, :], dims=2))

    num_t = size(data, 1)
    x_itp = LinearInterpolation(LinRange(0, 1, num_t - t0 + 1), [data[t, x_c:2:end-1, :] for t = t0:num_t])

    Q1_itp = LinearInterpolation(LinRange(0, 1, num_t - t0 + 1), [Q[t, 1, :, :] for t = t0:num_t])
    Qend_itp = LinearInterpolation(LinRange(0, 1, num_t - t0 + 1), [Q[t, end, :, :] for t = t0:num_t])

    # boundary polyhedron definition
    normals = vars["n"]
    ds = vars["d"]

    # compute growth rate
    num_s = size(x, 1)
    segment_lengths = [sqrt.(sum(abs2, data[t, x_c:2:end-1, :][2:end, :] - data[t, x_c:2:end-1, :][1:end-1, :]; dims=2)) for t = t0:num_t]
    arclengths = [cumsum([0.0; segment_lengths[t]], dims=1) for t = 1:length(segment_lengths)]
    length_t = [arclengths[t][end] for t = 1:length(arclengths)]
    fractional = [arclengths[t] ./ length_t[t] for t = 1:length(arclengths)]
    # smooth the growth rate - this is also reducing the magnitude.
    for k = 1:4
        segment_lengths[2:end-1] += 0.5 * (segment_lengths[1:end-2] + segment_lengths[3:end] - 2 * segment_lengths[2:end-1])
    end
    delta_lengths = [segment_lengths[t+1] - segment_lengths[t] for t = 1:length(segment_lengths)-1]
    delta_lengths_scoord = [LinearInterpolation(vec(fractional[t][1:end-1]), vec(delta_lengths[t]); extrapolation_bc=Flat()) for t = 1:length(delta_lengths)]

    growth_arr = zeros(num_t - t0, num_s)
    s_i = [s for s = LinRange(0, 1, num_s)]
    for t = 1:(num_t-t0)
        growth_arr[t, :] = delta_lengths_scoord[t](s_i)
    end

    dt_growth = 1.0 / (num_t - t0)

    ## interpolate circumradii
    c_radii = vars["circum_radii"][t0:end, x_c:2:end-1]
    cradii_arr = zeros(num_t - t0, num_s)
    c_radii_scoord = [LinearInterpolation(vec(fractional[t][1:end]), vec(c_radii[t, 1:end]); extrapolation_bc=Flat()) for t = 1:length(delta_lengths)]
    for t = 1:(num_t-t0)
        cradii_arr[t, :] = 1.0 * c_radii_scoord[t](s_i)
    end

    return (; sample_id, cfg, x, Qstart, x_itp, Q1_itp, Qend_itp,
            normals, ds, num_t, num_s, t0, growth_arr, cradii_arr, dt_growth, s_i)
end


# ---- per-run parameters ----------------------------------------------------------------
function prepare_run(smp, T::Float64, twist_ratio::Float64, plastic_rate::Float64)
    x = copy(smp.x)                 # local copy -- safe under parallelism
    asym = 1.0
    x[:, 2] = asym * x[:, 2]
    max_x = maximum(x[:, 1])

    N = size(x, 1) - 1
    Ns = N

    L = sum([norm(x[i+1, :] - x[i, :]) for i = 1:size(x, 1)-1])

    x0 = copy(x)
    x0 = x0 / L

    dxmax = 2.0 * maximum([norm(x0[i+1, :] - x0[i, :]) for i = 1:size(x, 1)-1])

    walls = make_walls(smp.normals, smp.ds / L) # need to rescale d by length unit

    S = 1e3
    St = 1e3
    B = 1.0
    Bt = twist_ratio * B
    alpha = 3e2 # self-contact parameter

    gamma_p = 1.0
    gamma_t = 1.0

    dx = minimum([norm(x0[i+1, :] - x0[i, :]) for i = 1:size(x0, 1)-1])
    x_start = copy(x0)
    Qstart = copy(smp.Qstart)       # local copy -- run_simulation normalizes in place

    dt = 3e-4 * dx
    n_saves = 200

    ## make growth rate interpolant
    # the measured interval is mapped onto [0, T], so the growth rate is rescaled
    # by T: the total accumulated growth is the same for every T in the sweep.
    nodes = ([t for t = LinRange(0, T, smp.num_t - smp.t0)], [s for s = LinRange(0, 1, smp.num_s)])
    growth_rate = extrapolate(interpolate(nodes, smp.growth_arr / (L * smp.dt_growth * T), Gridded(Linear())), Flat())
    circum_radii = extrapolate(interpolate(nodes, smp.cradii_arr / L, Gridded(Linear())), Flat())

    W1 = 10000.0   # wall stiffness of ext_forces
    cutoff = 0.1
    params = [T, plastic_rate, S, St, alpha, twist_ratio]

    return (; x_start, Qstart, dt, Ns, N, L, T, n_saves, S, St, B, Bt, gamma_p, gamma_t, alpha,
              walls, cutoff, W1, growth_rate, circum_radii, plastic_rate, dxmax, params)
end

# ---- clamped, data-driven ends (fast layout, in place) -----------------------------------
function make_bc_fast(smp, T::Float64, L::Float64)
    x_itp = smp.x_itp
    Q1_itp = smp.Q1_itp
    Qend_itp = smp.Qend_itp
    function bc!(t, x, Q)
        Q[1] = M3(Q1_itp(t / T))
        Q[end] = M3(Qend_itp(t / T))
        xi = x_itp(t / T)
        x[1] = V3(xi[1, :] / L)
        x[end] = V3(xi[end, :] / L)
        return nothing
    end
    return bc!
end
