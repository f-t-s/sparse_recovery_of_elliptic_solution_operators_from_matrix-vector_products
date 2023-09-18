using LinearAlgebra
using MKLSparse
using CompressingSolvers
using DelimitedFiles
using Random
using ArgParse
include("coefficients_and_potentials.jl")

settings = ArgParseSettings()
@add_arg_table settings begin
    "--model_size"
        help = "Which of the three models. Small(1), Medium(2), or Large(3)."
        required = true
        arg_type = Int
end
parsed_args = parse_args(settings)

model_size = parsed_args["model_size"]

@show model_size

pb = gridap_poisson("./gridap_models/demo-$(model_size).json")


ρ_list = [3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
num_threads = Threads.nthreads()
BLAS.set_num_threads(num_threads)

# Just for compilation
t_geo_list = Float64[]
t_meas_list = Float64[]
n_meas_list = Int64[]
t_rec_list = Float64[]
rel_err_list = Float64[]
for ρ in [minimum(ρ_list)]
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
    @show ρ
    # deleting the former reconstruction 
    rk = 0
    GC.gc()
    rk, logs = reconstruct(pb, ρ)
    push!(t_geo_list, logs[:t_geo])
    push!(t_meas_list, logs[:t_meas])
    push!(n_meas_list, logs[:n_meas])
    push!(t_rec_list, logs[:t_rec])
    push!(rel_err_list, CompressingSolvers.compute_relative_error(rk, pb))
end

base_path = "./out/csv/gridap/n_threads_$(num_threads)_model_size_$(model_size)"
writedlm(base_path * "_t_geo.csv", t_geo_list, ',')
writedlm(base_path * "_t_meas.csv", t_meas_list, ',')
writedlm(base_path * "_n_meas.csv", n_meas_list, ',')
writedlm(base_path * "_t_rec.csv", t_rec_list, ',')
writedlm(base_path * "_rel_err.csv", rel_err_list, ',')
writedlm(base_path * "_rho_list.csv", ρ_list, ',')
