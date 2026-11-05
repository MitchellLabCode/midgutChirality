## Production runs: every sample × the parameter grid (T, twist_ratio, plastic_rate) on the V6 data
#
# Measured samples (prod_common.jl) with their walls, self-contact, remeshing, data-driven clamped
# ends, growth and remodeling, with
#   - refIV6_sample*_walls3_cball.mat as reference data (see prod_common.jl)
#   - the fast kernel (fastprod.jl + ../src/fastkernel.jl), checked against
#     run_simulation_plastic_remeshed (../src/rod_contact2.jl) by validate_fastprod.jl
#   - rest curvature in the material frame (material_u0 = true) and growth of remeshed edges
#     weighted by w_seg, as in ../src/rod_contact2.jl
#
# Usage (from this folder):
#   julia -p 48 QTpTw_batch_samples_V6.jl          # runs `default_samples`
#   julia -p 48 QTpTw_batch_samples_V6.jl 4 5 6    # runs only samples 4, 5 and 6
#
# Output: ./<date>_QTpTw_contact2_t<t0>_w3_a3e2_s<id>V6/cutoff01contact_T<T>_tw<tw>_p<p>.mat
# with keys xsaves, Qsaves, params. Runs are resumable.

using Distributed

@everywhere begin import Pkg; Pkg.activate("../src/rods")
    Pkg.instantiate();
end
@everywhere begin # imports
    include("../src/rod_contact2.jl")
    using LinearAlgebra
    using Random
    using MAT
    using Interpolations
    include("../src/fastkernel.jl")
    include("fastprod.jl")
    include("prod_common.jl")
end

@everywhere begin
    function run_one(smp, T::Float64, twist_ratio::Float64, plastic_rate::Float64)
        folder_name = "./$(date_str)_QTpTw_contact2_t$(smp.t0)_w3_a3e2_$(smp.cfg.tag)"
        out_filename = folder_name * "/cutoff01contact_T$(T)_tw$(twist_ratio)_p$(plastic_rate).mat"
        if isfile(out_filename)
            return out_filename
        end
        mkpath(folder_name)

        r = prepare_run(smp, T, twist_ratio, plastic_rate)
        bc! = make_bc_fast(smp, T, r.L)
        println("sample $(smp.sample_id): T=$(T) tw=$(twist_ratio) p=$(plastic_rate) dt=$(r.dt)")
        flush(stdout)

        t_start = time()
        xsaves, Qsaves = run_fastprod(r.x_start, r.Qstart, r.dt, r.Ns, T, r.n_saves,
            r.S, r.St, r.B, r.Bt, r.gamma_p, r.gamma_t, bc!, r.circum_radii, r.alpha, r.walls, r.cutoff,
            r.W1, r.growth_rate, plastic_rate, r.dxmax; material_u0=true)
        elapsed = time() - t_start

        matwrite(out_filename, Dict("xsaves" => xsaves, "Qsaves" => Qsaves, "params" => r.params,
                                    "ref" => smp.cfg.ref, "material_u0" => true, "elapsed_s" => elapsed))
        return out_filename
    end

    # per-worker cache of the sample setups
    sample_cache = Dict{Int,Any}()
    get_sample(id) = get!(() -> setup_sample(id), sample_cache, id)
end


# ---- sample selection: command-line ids override `default_samples` ---------------
samples_to_run = isempty(ARGS) ? default_samples : parse.(Int, ARGS)
for s in samples_to_run
    haskey(sample_configs, s) || error("Unknown sample id $(s). Known ids: $(sort(collect(keys(sample_configs))))")
end

# all (sample, T, tw, p) jobs in a single pmap, shuffled once on the driver for load
# balancing (the T = 3 runs take about 3x longer than T = 1)
jobs = shuffle([(s, v...) for s in samples_to_run for v in values_vec])

println("---- Grid: $(length(Ts)) T x $(length(twist_ratios)) twist x $(length(plastic_rates_arr)) plastic = $(length(values_vec)) runs/sample ----")
println("---- Samples: $(samples_to_run) -> $(length(values_vec) * length(samples_to_run)) simulations total ----")

status = pmap(jobs) do j
    sample_id, T, twist_ratio, plastic_rate = j
    run_one(get_sample(sample_id), T, twist_ratio, plastic_rate)
end
println("~success~")
