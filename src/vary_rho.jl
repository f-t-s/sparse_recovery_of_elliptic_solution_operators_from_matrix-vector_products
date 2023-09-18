using LinearAlgebra
using MKLSparse
using CompressingSolvers
using DelimitedFiles
using Random
const CS = CompressingSolvers

q_list = [4, 5, 6, 7, 8, 9]
ρ_list = [4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]

num_threads = Threads.nthreads()
BLAS.set_num_threads(num_threads)

# create one seed to be used throughout to ensure that for a given q,
# all different values of rho are applied to the same problem.
sd = rand(1 : 10^6)

# Just for compilation
for q in [5]
    t_list = Float64[]
    rel_err_list = Float64[]
    n_matvec_list = Int[]
    for ρ in [5.0]
        @show ρ
        # resetting random seed
        Random.seed!(sd)
        A, coarse_domains, scales, basis_functions, multicolor_ordering, fine_domains, tree_function =
            CS.FD_periodic_Laplacian_subdivision_2d(q, ρ)
        measurement_matrix = CS.form_measurement_matrix(multicolor_ordering)
        measurement_results = CS.measure(cholesky(A), measurement_matrix)
        @time t = @elapsed L = CS.reconstruct(multicolor_ordering, 
                                        CS.center.(fine_domains), 
                                        measurement_matrix,
                                        measurement_results,
                                        tree_function) 
        @time rel_err = CompressingSolvers.compute_relative_error(L, cholesky(A), 200, 10)
        n_matvec = size(measurement_matrix, 2)
        push!(n_matvec_list, n_matvec)
        push!(rel_err_list, rel_err)
        push!(t_list, t)
    end
end

for q in q_list
    t_list = Float64[]
    rel_err_list = Float64[]
    n_matvec_list = Int[]
    for ρ in ρ_list
        @show ρ
        # resetting random seed
        Random.seed!(sd)
        A, coarse_domains, scales, basis_functions, multicolor_ordering, fine_domains, tree_function =
            CS.FD_periodic_Laplacian_subdivision_2d(q, ρ)
        measurement_matrix = CS.form_measurement_matrix(multicolor_ordering)
        measurement_results = CS.measure(cholesky(A), measurement_matrix)
        @time t = @elapsed L = CS.reconstruct(multicolor_ordering, 
                                        CS.center.(fine_domains), 
                                        measurement_matrix,
                                        measurement_results,
                                        tree_function) 
        @time rel_err = CompressingSolvers.compute_relative_error(L, cholesky(A), 200, 10)
        n_matvec = size(measurement_matrix, 2)
        push!(n_matvec_list, n_matvec)
        push!(rel_err_list, rel_err)
        push!(t_list, t)
    end
    writedlm("./out/csv/n_threads_$(num_threads)_matvecs_q_$(q)_vary_rho.csv", n_matvec_list, ',')
    writedlm("./out/csv/n_threads_$(num_threads)_rel_err_q_$(q)_vary_rho.csv", rel_err_list, ',')
    writedlm("./out/csv/n_threads_$(num_threads)_t_q_$(q)_vary_rho.csv", t_list, ',')
    writedlm("./out/csv/n_threads_$(num_threads)_rho_list_q_$(q).csv",ρ_list, ',')
    writedlm("./out/csv/matvecs_vs_rho_q_$(q).csv", hcat(ρ_list, n_matvec_list), ' ')
    writedlm("./out/csv/err_vs_rho_q_$(q).csv", hcat(ρ_list, rel_err_list), ' ')
end