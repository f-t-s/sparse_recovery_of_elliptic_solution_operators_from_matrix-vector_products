using LinearAlgebra
using MKLSparse
using CompressingSolvers
using DelimitedFiles
using Random
using ArgParse
include("coefficients_and_potentials.jl")

settings = ArgParseSettings()
@add_arg_table settings begin
    "--n"
        help = "The number of degrees of freedom per dimension"
        required = true
        arg_type = Int
    "--d" 
        help = "The spatial dimension"
        required = true
        arg_type = Int
    "--s"
        help = "Chooses the (fractional) order of the problem"
        required = true 
        arg_type = Float64
end
parsed_args = parse_args(settings)

n = parsed_args["n"]
d = parsed_args["d"]
s = parsed_args["s"]

@show n
@show d
@show s

# create the problem
if d == 2 
    ρ_list = [4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    pb = uniform2d_fractional(n, s, 1.0)
elseif d == 3
    ρ_list = [3.0, 4.0, 5.0, 6.0]
    pb = uniform3d_fractional(n, s, 1.0)
end
num_threads = Threads.nthreads()
BLAS.set_num_threads(num_threads)

# Just for compilation
t_geo_list = Float64[]
t_meas_list = Float64[]
n_meas_list = Int64[]
t_rec_list = Float64[]
rel_err_list = Float64[]
println("Running a test run in order to precompile")
for ρ in [minimum(ρ_list)]
    GC.gc()
    rk, logs = reconstruct(pb, ρ)
    push!(t_geo_list, logs[:t_geo])
    push!(t_meas_list, logs[:t_meas])
    push!(n_meas_list, logs[:n_meas])
    push!(t_rec_list, logs[:t_rec])
    push!(rel_err_list, CompressingSolvers.compute_relative_error(rk, pb))
end

# The actual computation
t_geo_list = Float64[]
t_meas_list = Float64[]
n_meas_list = Int64[]
t_rec_list = Float64[]
rel_err_list = Float64[]
for ρ in ρ_list
    println("Beginning run with rho=$(ρ)")
    # unallocating the previous solution to save memory
    rk = 0
    GC.gc()
    @show ρ
    rk, logs = reconstruct(pb, ρ)
    push!(t_geo_list, logs[:t_geo])
    push!(t_meas_list, logs[:t_meas])
    push!(n_meas_list, logs[:n_meas])
    push!(t_rec_list, logs[:t_rec])
    push!(rel_err_list, CompressingSolvers.compute_relative_error(rk, pb))
end

base_path = "./out/csv/fractional_poisson/n_threads_$(num_threads)_n_$(n)_d_$(d)_s_$(s)"
writedlm(base_path * "_t_geo.csv", t_geo_list, ',')
writedlm(base_path * "_t_meas.csv", t_meas_list, ',')
writedlm(base_path * "_n_meas.csv", n_meas_list, ',')
writedlm(base_path * "_t_rec.csv", t_rec_list, ',')
writedlm(base_path * "_rel_err.csv", rel_err_list, ',')
writedlm(base_path * "_rho_list.csv", ρ_list, ',')
