using DelimitedFiles
# This file compiles the results of individual runs 

# for n in [64; 128; 256; 512; 1024]
for model_size in 1 : 2
    path = "./out/csv/gridap/n_threads_24_model_size_$(model_size)_"
    ρs = readdlm(path * "rho_list.csv", ',')[:]
    n_meas = readdlm(path * "n_meas.csv", ',')[:]
    rel_errs = readdlm(path * "rel_err.csv", ',')[:]
    writedlm("./figures/gridap/$(model_size).csv", hcat(ρs, n_meas, rel_errs), ',')
end