using DelimitedFiles
# This file compiles the results of individual runs 

# Dirichlet Laplacian 

boundary="dirichlet"
potential="random_mild"
coefficients="constant"
d = "2"

for n in [64; 128; 256; 512; 1024]
    path = "./out/csv/fd_poisson/n_threads_24_n_$(n)_d_$(d)_boundary_$(boundary)_coefficients_$(coefficients)_potential_$(potential)_"
    ρs = readdlm(path * "rho_list.csv", ',')[:]
    n_meas = readdlm(path * "n_meas.csv", ',')[:]
    rel_errs = readdlm(path * "rel_err.csv", ',')[:]
    writedlm("./figures/$(n)_$(d)_$(boundary)_$(coefficients)_$(potential).csv", hcat(ρs, n_meas, rel_errs), ',')
end

boundary="periodic"
potential="random_mild"
coefficients="constant"
d = "2"

for n in [64; 128; 256; 512; 1024]
    path = "./out/csv/fd_poisson/n_threads_24_n_$(n)_d_$(d)_boundary_$(boundary)_coefficients_$(coefficients)_potential_$(potential)_"
    ρs = readdlm(path * "rho_list.csv", ',')[:]
    n_meas = readdlm(path * "n_meas.csv", ',')[:]
    rel_errs = readdlm(path * "rel_err.csv", ',')[:]
    writedlm("./figures/$(n)_$(d)_$(boundary)_$(coefficients)_$(potential).csv", hcat(ρs, n_meas, rel_errs), ',')
end

boundary="neumann"
potential="random_mild"
coefficients="constant"
d = "2"

for n in [64; 128; 256; 512; 1024]
    path = "./out/csv/fd_poisson/n_threads_24_n_$(n)_d_$(d)_boundary_$(boundary)_coefficients_$(coefficients)_potential_$(potential)_"
    ρs = readdlm(path * "rho_list.csv", ',')[:]
    n_meas = readdlm(path * "n_meas.csv", ',')[:]
    rel_errs = readdlm(path * "rel_err.csv", ',')[:]
    writedlm("./figures/$(n)_$(d)_$(boundary)_$(coefficients)_$(potential).csv", hcat(ρs, n_meas, rel_errs), ',')
end

boundary="dirichlet"
potential="random_mild"
coefficients="random_severe"
d = "2"

for n in [64; 128; 256; 512; 1024]
    path = "./out/csv/fd_poisson/n_threads_24_n_$(n)_d_$(d)_boundary_$(boundary)_coefficients_$(coefficients)_potential_$(potential)_"
    ρs = readdlm(path * "rho_list.csv", ',')[:]
    n_meas = readdlm(path * "n_meas.csv", ',')[:]
    rel_errs = readdlm(path * "rel_err.csv", ',')[:]
    writedlm("./figures/$(n)_$(d)_$(boundary)_$(coefficients)_$(potential).csv", hcat(ρs, n_meas, rel_errs), ',')
end

boundary="periodic"
potential="random_mild"
coefficients="random_severe"
d = "2"

for n in [64; 128; 256; 512; 1024]
    path = "./out/csv/fd_poisson/n_threads_24_n_$(n)_d_$(d)_boundary_$(boundary)_coefficients_$(coefficients)_potential_$(potential)_"
    ρs = readdlm(path * "rho_list.csv", ',')[:]
    n_meas = readdlm(path * "n_meas.csv", ',')[:]
    rel_errs = readdlm(path * "rel_err.csv", ',')[:]
    writedlm("./figures/$(n)_$(d)_$(boundary)_$(coefficients)_$(potential).csv", hcat(ρs, n_meas, rel_errs), ',')
end

boundary="neumann"
potential="random_mild"
coefficients="random_severe"
d = "2"

for n in [64; 128; 256; 512; 1024]
    path = "./out/csv/fd_poisson/n_threads_24_n_$(n)_d_$(d)_boundary_$(boundary)_coefficients_$(coefficients)_potential_$(potential)_"
    ρs = readdlm(path * "rho_list.csv", ',')[:]
    n_meas = readdlm(path * "n_meas.csv", ',')[:]
    rel_errs = readdlm(path * "rel_err.csv", ',')[:]
    writedlm("./figures/$(n)_$(d)_$(boundary)_$(coefficients)_$(potential).csv", hcat(ρs, n_meas, rel_errs), ',')
end

################################
# d = 3
################################


boundary="dirichlet"
potential="random_mild"
coefficients="constant"
d = "3"

for n in [32; 64]
    path = "./out/csv/fd_poisson/n_threads_24_n_$(n)_d_$(d)_boundary_$(boundary)_coefficients_$(coefficients)_potential_$(potential)_"
    ρs = readdlm(path * "rho_list.csv", ',')[:]
    n_meas = readdlm(path * "n_meas.csv", ',')[:]
    rel_errs = readdlm(path * "rel_err.csv", ',')[:]
    writedlm("./figures/$(n)_$(d)_$(boundary)_$(coefficients)_$(potential).csv", hcat(ρs, n_meas, rel_errs), ',')
end

boundary="periodic"
potential="random_mild"
coefficients="constant"
d = "3"

for n in [32; 64]
    path = "./out/csv/fd_poisson/n_threads_24_n_$(n)_d_$(d)_boundary_$(boundary)_coefficients_$(coefficients)_potential_$(potential)_"
    ρs = readdlm(path * "rho_list.csv", ',')[:]
    n_meas = readdlm(path * "n_meas.csv", ',')[:]
    rel_errs = readdlm(path * "rel_err.csv", ',')[:]
    writedlm("./figures/$(n)_$(d)_$(boundary)_$(coefficients)_$(potential).csv", hcat(ρs, n_meas, rel_errs), ',')
end

boundary="neumann"
potential="random_mild"
coefficients="constant"
d = "3"

for n in [32; 64]
    path = "./out/csv/fd_poisson/n_threads_24_n_$(n)_d_$(d)_boundary_$(boundary)_coefficients_$(coefficients)_potential_$(potential)_"
    ρs = readdlm(path * "rho_list.csv", ',')[:]
    n_meas = readdlm(path * "n_meas.csv", ',')[:]
    rel_errs = readdlm(path * "rel_err.csv", ',')[:]
    writedlm("./figures/$(n)_$(d)_$(boundary)_$(coefficients)_$(potential).csv", hcat(ρs, n_meas, rel_errs), ',')
end

boundary="dirichlet"
potential="random_mild"
coefficients="random_severe"
d = "3"

for n in [32; 64]
    path = "./out/csv/fd_poisson/n_threads_24_n_$(n)_d_$(d)_boundary_$(boundary)_coefficients_$(coefficients)_potential_$(potential)_"
    ρs = readdlm(path * "rho_list.csv", ',')[:]
    n_meas = readdlm(path * "n_meas.csv", ',')[:]
    rel_errs = readdlm(path * "rel_err.csv", ',')[:]
    writedlm("./figures/$(n)_$(d)_$(boundary)_$(coefficients)_$(potential).csv", hcat(ρs, n_meas, rel_errs), ',')
end

boundary="periodic"
potential="random_mild"
coefficients="random_severe"
d = "3"

for n in [32; 64]
    path = "./out/csv/fd_poisson/n_threads_24_n_$(n)_d_$(d)_boundary_$(boundary)_coefficients_$(coefficients)_potential_$(potential)_"
    ρs = readdlm(path * "rho_list.csv", ',')[:]
    n_meas = readdlm(path * "n_meas.csv", ',')[:]
    rel_errs = readdlm(path * "rel_err.csv", ',')[:]
    writedlm("./figures/$(n)_$(d)_$(boundary)_$(coefficients)_$(potential).csv", hcat(ρs, n_meas, rel_errs), ',')
end

boundary="neumann"
potential="random_mild"
coefficients="random_severe"
d = "3"

for n in [32; 64]
    path = "./out/csv/fd_poisson/n_threads_24_n_$(n)_d_$(d)_boundary_$(boundary)_coefficients_$(coefficients)_potential_$(potential)_"
    ρs = readdlm(path * "rho_list.csv", ',')[:]
    n_meas = readdlm(path * "n_meas.csv", ',')[:]
    rel_errs = readdlm(path * "rel_err.csv", ',')[:]
    writedlm("./figures/$(n)_$(d)_$(boundary)_$(coefficients)_$(potential).csv", hcat(ρs, n_meas, rel_errs), ',')
end