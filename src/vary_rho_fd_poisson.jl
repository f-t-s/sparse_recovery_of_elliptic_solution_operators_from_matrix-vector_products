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
    "--boundary"
        help = "Chooses the boundary conditions"
        required = true 
        arg_type = String
    "--coefficients"
        help = "Chooses the conductivity coefficients"
        required = true 
        arg_type = String
    "--potential"
        help = "Chooses the potential"
        required = true 
        arg_type = String
end
parsed_args = parse_args(settings)

n = parsed_args["n"]
d = parsed_args["d"]
boundary = parsed_args["boundary"]
coefficients = parsed_args["coefficients"] 
potential = parsed_args["potential"]

@show n
@show d
@show boundary
@show coefficients
@show potential

α, β = select_coefficients_and_potential(d, coefficients, potential)

# create the problem
if d == 2 
    ρ_list = [4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    if boundary == "periodic"
        pb = uniform2d_periodic_fd_poisson(n, α, β)
    elseif boundary == "dirichlet"
        pb = uniform2d_dirichlet_fd_poisson(n, α, β)
    elseif boundary == "neumann"
        pb = uniform2d_neumann_fd_poisson(n, α, β)
    end
elseif d == 3
    ρ_list = [3.0, 4.0, 5.0, 6.0]
    if boundary == "periodic"
        pb = uniform3d_periodic_fd_poisson(n, α, β)
    elseif boundary == "dirichlet"
        pb = uniform3d_dirichlet_fd_poisson(n, α, β)
    elseif boundary == "neumann"
        pb = uniform3d_neumann_fd_poisson(n, α, β)
    end
end
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

base_path = "./out/csv/fd_poisson/n_threads_$(num_threads)_n_$(n)_d_$(d)_boundary_$(boundary)_coefficients_$(coefficients)_potential_$(potential)"
writedlm(base_path * "_t_geo.csv", t_geo_list, ',')
writedlm(base_path * "_t_meas.csv", t_meas_list, ',')
writedlm(base_path * "_n_meas.csv", n_meas_list, ',')
writedlm(base_path * "_t_rec.csv", t_rec_list, ',')
writedlm(base_path * "_rel_err.csv", rel_err_list, ',')
writedlm(base_path * "_rho_list.csv", ρ_list, ',')
