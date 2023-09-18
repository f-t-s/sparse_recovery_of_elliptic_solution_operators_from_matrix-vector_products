using DelimitedFiles
# This file compiles the results of individual runs 

s=0.5
d = "2"

# for n in [64; 128; 256; 512; 1024]
for n in [64; 128; 256; 512]
    path = "./out/csv/fractional_poisson/n_threads_24_n_$(n)_d_$(d)_s_$(s)_"
    ρs = readdlm(path * "rho_list.csv", ',')[:]
    n_meas = readdlm(path * "n_meas.csv", ',')[:]
    rel_errs = readdlm(path * "rel_err.csv", ',')[:]
    writedlm("./figures/fractional/$(n)_$(d)_$(s).csv", hcat(ρs, n_meas, rel_errs), ',')
end

s=1.5
d = "2"

# for n in [64; 128; 256; 512; 1024]
for n in [64; 128; 256; 512]
    path = "./out/csv/fractional_poisson/n_threads_24_n_$(n)_d_$(d)_s_$(s)_"
    ρs = readdlm(path * "rho_list.csv", ',')[:]
    n_meas = readdlm(path * "n_meas.csv", ',')[:]
    rel_errs = readdlm(path * "rel_err.csv", ',')[:]
    writedlm("./figures/fractional/$(n)_$(d)_$(s).csv", hcat(ρs, n_meas, rel_errs), ',')
end

s=0.5
d = "3"

for n in [32; 64;]
    path = "./out/csv/fractional_poisson/n_threads_24_n_$(n)_d_$(d)_s_$(s)_"
    ρs = readdlm(path * "rho_list.csv", ',')[:]
    n_meas = readdlm(path * "n_meas.csv", ',')[:]
    rel_errs = readdlm(path * "rel_err.csv", ',')[:]
    writedlm("./figures/fractional/$(n)_$(d)_$(s).csv", hcat(ρs, n_meas, rel_errs), ',')
end

s=1.5
d = "3"

for n in [32; 64;]
    path = "./out/csv/fractional_poisson/n_threads_24_n_$(n)_d_$(d)_s_$(s)_"
    ρs = readdlm(path * "rho_list.csv", ',')[:]
    n_meas = readdlm(path * "n_meas.csv", ',')[:]
    rel_errs = readdlm(path * "rel_err.csv", ',')[:]
    writedlm("./figures/fractional/$(n)_$(d)_$(s).csv", hcat(ρs, n_meas, rel_errs), ',')
end
