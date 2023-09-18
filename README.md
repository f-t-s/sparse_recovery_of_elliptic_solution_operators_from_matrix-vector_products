# sparse_recovery_of_elliptic_solution_operators_from_matrix-vector_products
This repository contains the code for reproducing the experiments in the paper 
[Sparse recovery of elliptic solution operators from matrix-vector products](https://arxiv.org/pdf/2110.05351.pdf),
to appear in [SISC](https://www.siam.org/publications/journals/siam-journal-on-scientific-computing-sisc). 

They rely on the (publicly available) custom Julia package [CompressingSolvers.jl](https://github.com/f-t-s/CompressingSolvers.jl). The environment
is specified using the files ``Manifest.toml`` and ``Project.toml``.

The experiments in Figures 7 and 8 were run using the script ``fd_poisson.pbs``, that calls 
``julia --threads=24 --project=. ./src/vary_rho_fd_poisson.jl --n $N --d $D --boundary $BOUNDARY --coefficients $COEFFICIENTS --potential $POTENTIAL``.
The parameter ``d`` is the spatial dimension. ``n`` is the number of grid points in each dimension, resulting in a total of n^d grid points. 
The parameters ``boundary``, ``coefficients``, ``potential`` prescribe the boundary conditions, the conductivity coefficient field, and the zero order,
potential term of the elliptic PDE. The plots in Figure 7 use ``boundary=periodic``, ``coefficients=constant``, and ``potential=random_mild``. Those
in Figure 8 use ``boundary=periodic``, ``coefficients=periodic_severe``, and ``potential=random_mild``.

The experiments in Figure 9 were run using the script ``fractional_poisson.pbs``, which calls 
``julia --threads=24 --project=. ./src/vary_rho_fractional_poisson.jl --n $N --d $D --s $S``. 
Here, ``n`` and ``d`` have the same meaning as before, and ``s`` specifies the order of the PDE. 
For instance, the Laplace operator corresponds to ``s=1``. 

Finally, the results in Figure 11 were run using ``gridap_poisson.pbs``, which calls 
``julia --threads=24 --project=. ./src/vary_rho_gridap_poisson.jl --model_size $MODEL_SIZE``. Here, ``model_size`` is an integer between 1 and 3 
that specifies the size of the FEM model to be used. 

