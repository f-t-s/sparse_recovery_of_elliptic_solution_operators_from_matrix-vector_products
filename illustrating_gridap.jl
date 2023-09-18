using Random
Random.seed!(172)
using LinearAlgebra: Matrix
using CompressingSolvers
using SparseArrays
using LinearAlgebra
using GLMakie
using GridapMakie
using Gridap
using FlorianStyle

# Setting up the test domains
ρ = 3
# ρ = Inf

# n = 2^7
# pb = uniform2d_dirichlet_fd_poisson(n)
# pb = uniform3d_dirichlet_fd_poisson(n)

# pb = uniform2d_neumann_fd_poisson(n)

# pb = uniform2d_fractional(n, 0.50, 1.0)
# pb = uniform2d_periodic_fd_poisson(n)
# pb = uniform2d_dirichlet_fd_poisson(n)

model_size=1

path1 = "./gridap_models/demo-1.json"
path2 = "./gridap_models/demo-2.json"
path3 = "./gridap_models/demo-3.json"

# path = "./gridap_models/demo_refined.msh"
# # path = "./gridap_models/model.json"

model = DiscreteModelFromFile(path1)
Ω = Triangulation(model)
∂Ω = BoundaryTriangulation(model)
Γ = BoundaryTriangulation(model, tags=["boundary1", "boundary2"])
fig = plot(Ω, shading=true, color=silver, axis=(show_axis=false,), figure=(resolution=(1800,1200),))
wireframe!(∂Ω, color=darksky)
wireframe!(Γ, color=joshua)
# save("./out/figures/gridap_full_mesh.png", fig)

## model1 = DiscreteModelFromFile(path1)
## ∂Ω1 = BoundaryTriangulation(model1)
## Γ1 = BoundaryTriangulation(model1, tags=["boundary1", "boundary2"])
## 
## 
## model2 = DiscreteModelFromFile(path2)
## # fig = wireframe(Ω, linewidth=0.5, shading=true, color=silver, axis=(show_axis=false,), figure=(resolution=(1800,1200),))
## ∂Ω2 = BoundaryTriangulation(model2)
## Γ2 = BoundaryTriangulation(model2, tags=["boundary1", "boundary2"])
## 
## model3 = DiscreteModelFromFile(path3)
## ∂Ω3 = BoundaryTriangulation(model3)
## Γ3 = BoundaryTriangulation(model3, tags=["boundary1", "boundary2"])
## 
## Ω = Triangulation(model3)
## 
## # fig = plot(Ω, shading=true, color=silver, axis=(show_axis=false,), figure=(resolution=(1800,1200),))
## fig = plot(Ω, shading=true, color=silver, axis=(show_axis=false,), figure=(resolution=(1800,1200),))
## 
## wireframe!(∂Ω3, linewidth=2.0, color=steelblue)
## # wireframe!(∂Ω3, linewidth=2.0, color=darksky, overdraw=true)
## # wireframe!(∂Ω2, linewidth=4.0, color=seagreen, overdraw=true)
## # wireframe!(∂Ω1, linewidth=8.0, color=steelblue, overdraw=true)
## wireframe!(Γ3, linewidth=2.0, color=rust, overdraw=true)
## wireframe!(Γ2, linewidth=4.0, color=orange, overdraw=true)
## wireframe!(Γ1, linewidth=8.0, color=joshua, overdraw=true)
## 
## plot!(Ω, color=silver)
# wireframe!(Ω, color =)

# save("./out/figures/gridap_partial_mesh.png", fig)




# save("./out/figures/gridap_zoomed_mesh_$(model_size).png", fig)

