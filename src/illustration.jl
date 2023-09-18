using LinearAlgebra
using CairoMakie
using DelimitedFiles
using FlorianStyle
using Colors
set_theme!(empty_theme)

q = 6
n = 2 ^ q
N = n ^ 2 
h = 1 / (n + 1)

levels = [(4^k + 1) : 4^(k + 1) for k = 1 : (q - 1)]
levels = vcat([1 : 4], levels)

li = LinearIndices((n, n))

x = zeros(2, N)
A = zeros(N, N)  

for i = 1 : n, j = 1 : n
    x[1, li[i, j]] = i * h
    x[2, li[i, j]] = j * h

    A[li[i,j], li[i,j]] = 4 / h 
    if i > 1
        A[li[i, j], li[i - 1, j]] = - 1/h
    end
    if i < n
        A[li[i, j], li[i + 1, j]] = - 1/h
    end

    if j > 1 A[li[i, j], li[i, j - 1]] = - 1/h
    end
    if j < n
        A[li[i, j], li[i, j + 1]] = - 1/h
    end
end

ws = Vector{Float64}[]
push!(ws, ones(N))
for k = 1 : q 
    bs = 2^(q - k) 
    for ic = 1 : 2^(k - 1), jc = 1 : 2^( k - 1) 
        push!(ws, zeros(N))
        last(ws)[li[((ic - 1) * 2 * bs + 1) : (ic - 1) * 2 * bs + bs, ((jc - 1) * 2 * bs + 1) : (jc - 1) * 2 * bs + bs]] .= 1.0
        last(ws)[li[((ic - 1) * 2 * bs + bs + 1) : (ic - 1) * 2 * bs + 2 * bs, ((jc - 1) * 2 * bs + 1) : (jc - 1) * 2 * bs + bs]] .= 1.0
        last(ws)[li[((ic - 1) * 2 * bs + 1) : (ic - 1) * 2 * bs + bs, ((jc - 1) * 2 * bs + bs + 1) : (jc - 1) * 2 * bs + 2 * bs]] .= - 1.0
        last(ws)[li[((ic - 1) * 2 * bs + bs + 1) : (ic - 1) * 2 * bs + 2 * bs, ((jc - 1) * 2 * bs + bs + 1) : (jc - 1) * 2 * bs + 2 * bs]] .= - 1.0

        push!(ws, zeros(N))
        last(ws)[li[((ic - 1) * 2 * bs + 1) : (ic - 1) * 2 * bs + bs, ((jc - 1) * 2 * bs + 1) : (jc - 1) * 2 * bs + bs]] .= 1.0
        last(ws)[li[((ic - 1) * 2 * bs + bs + 1) : (ic - 1) * 2 * bs + 2 * bs, ((jc - 1) * 2 * bs + 1) : (jc - 1) * 2 * bs + bs]] .= - 1.0
        last(ws)[li[((ic - 1) * 2 * bs + 1) : (ic - 1) * 2 * bs + bs, ((jc - 1) * 2 * bs + bs + 1) : (jc - 1) * 2 * bs + 2 * bs]] .= 1.0
        last(ws)[li[((ic - 1) * 2 * bs + bs + 1) : (ic - 1) * 2 * bs + 2 * bs, ((jc - 1) * 2 * bs + bs + 1) : (jc - 1) * 2 * bs + 2 * bs]] .= - 1.0

        push!(ws, zeros(N))
        last(ws)[li[((ic - 1) * 2 * bs + 1) : (ic - 1) * 2 * bs + bs, ((jc - 1) * 2 * bs + 1) : (jc - 1) * 2 * bs + bs]] .= 1.0
        last(ws)[li[((ic - 1) * 2 * bs + bs + 1) : (ic - 1) * 2 * bs + 2 * bs, ((jc - 1) * 2 * bs + 1) : (jc - 1) * 2 * bs + bs]] .= - 1.0
        last(ws)[li[((ic - 1) * 2 * bs + 1) : (ic - 1) * 2 * bs + bs, ((jc - 1) * 2 * bs + bs + 1) : (jc - 1) * 2 * bs + 2 * bs]] .= - 1.0
        last(ws)[li[((ic - 1) * 2 * bs + bs + 1) : (ic - 1) * 2 * bs + 2 * bs, ((jc - 1) * 2 * bs + bs + 1) : (jc - 1) * 2 * bs + 2 * bs]] .= 1.0
    end
end

W = reduce(hcat, ws)
W = W / diagm(sqrt.(diag(W' * W)))

L = Matrix(cholesky(Hermitian(inv(W' * A * W))).L)
WL = W * L

# surface(reshape(x[1, :], n, n), reshape(x[2, :], n, n), reshape(WL[:, 21], n, n))

# base_path = "./out/figures/decay_figure_q_$(q)"
# writedlm(base_path * "_L.csv", L)
# writedlm(base_path * "_WL.csv", WL)
# writedlm(base_path * "_W.csv", W)
# 
# writedlm(base_path * "_zero_column.csv", hcat(vec(x[1,:]), vec(x[2,:]), zeros(N)))
# for col in 1 : N
#     writedlm(base_path * "_WL_column_$(col).csv", hcat(vec(x[1,:]), vec(x[2,:]), WL[:, col] / maximum(abs.(WL[:, col]))))
#     writedlm(base_path * "_W_column_$(col).csv", hcat(vec(x[1,:]), vec(x[2,:]), W[:, col]) / maximum(abs.(W[:, col])))
# end
# 
# ct = CartesianIndices((N, N))
# out_L = Matrix(transpose(reduce(hcat, [[ct[k][1], ct[k][2], val] for (k, val) in enumerate(L)])))

# writedlm(base_path * "_coordinate_L.csv",out_L )
# writedlm(base_path * "_coordinate_log10_L.csv",hcat(out_L[:,1:2], log10.(1e-10 .+ abs.(out_L[:, end]))))


Lnorm = L ./ diag(L)

xlevels = Vector{Float64}[]
ylevels = Vector{Float64}[]
vlevels = Vector{Float64}[]
for k = 1 : q
    push!(xlevels, zeros(0))
    push!(ylevels, zeros(0))
    push!(vlevels, zeros(0))
    cart_inds = CartesianIndices((N, length(levels[k])))
    for (ind, val) in enumerate(Lnorm[:, levels[k]])
        push!(last(xlevels), cart_inds[ind][2] + (first(levels[k]) - 1))
        push!(last(ylevels), cart_inds[ind][1])
        push!(last(vlevels), val)
    end
end

if q == 6 
    fig = Figure()
    fig[1, 1] = ax = Axis(fig)
    ax.yreversed=true
    ax.aspect=AxisAspect(1)
    for k = 1 : q
        if k == 3 
            cmp = ["white", rust]
            heatmap!(ax, xlevels[k], ylevels[k], log10.(1e-10 .+ abs.(vlevels[k])), colormap=cmp)
        elseif k == 4
            cmp = ["white", steelblue]
            heatmap!(ax, xlevels[k], ylevels[k], log10.(1e-10 .+ abs.(vlevels[k])), colormap=cmp)
        elseif k == 5 
            cmp = ["white", joshua]
            heatmap!(ax, xlevels[k], ylevels[k], log10.(1e-10 .+ abs.(vlevels[k])), colormap=cmp)
        else
            cmp = ["white", seagreen]
            heatmap!(ax, xlevels[k], ylevels[k], log10.(1e-10 .+ abs.(vlevels[k])), colormap=cmp)
        end
    end

    for k = 1 : q
        writedlm("./out/figures/decay_figure/decay_figure_log10L_q_$(q)_k_$(k).csv", hcat(xlevels[k], ylevels[k], vlevels[k]))
    end
    save("./out/figures/decay_figure/decay_figure_q_$(q)_log10L.pdf", fig, resolution = (1000, 1000))
end


# Plotting a basis function on lvl 3
ind3 = 19
linds = LinearIndices((n, n))
xinds = 1 : Int(n / 2)
yinds = 1 : Int(n / 2)
# xinds = 1 : n
# yinds = 1 : n
xs = [x[1, linds[xind, yind]] for xind in xinds, yind in yinds][:]
ys = [x[2, linds[xind, yind]] for xind in xinds, yind in yinds][:]
vs = [WL[linds[xind, yind], ind3] for xind in xinds, yind in yinds][:]
ws = [W[linds[xind, yind], ind3] for xind in xinds, yind in yinds][:]
writedlm("./out/figures/decay_figure/adapted_basis_function_lvl_3_q_$(q).csv", hcat(xs, ys, vs))
writedlm("./out/figures/decay_figure/basis_function_lvl_3_q_$(q).csv", hcat(xs, ys, ws))


# Plotting a basis function on lvl 4
ind4 = 109
linds = LinearIndices((n, n))
xinds = 1 : Int(n / 2)
yinds = (Int(n / 2)) : n
# xinds = 1 : n
# yinds = 1 : n
xs = [x[1, linds[xind, yind]] for xind in xinds, yind in yinds][:]
ys = [x[2, linds[xind, yind]] for xind in xinds, yind in yinds][:]
vs = [WL[linds[xind, yind], ind4] for xind in xinds, yind in yinds][:]
ws = [W[linds[xind, yind], ind4] for xind in xinds, yind in yinds][:]
writedlm("./out/figures/decay_figure/adapted_basis_function_lvl_4_q_$(q).csv", hcat(xs, ys, vs))
writedlm("./out/figures/decay_figure/basis_function_lvl_4_q_$(q).csv", hcat(xs, ys, ws))

# Plotting a basis function on lvl 5
ind5 = 940
linds = LinearIndices((n, n))
xinds = (Int(n / 2) ) : n
yinds = 1 : Int(n / 2)
xs = [x[1, linds[xind, yind]] for xind in xinds, yind in yinds][:]
ys = [x[2, linds[xind, yind]] for xind in xinds, yind in yinds][:]
vs = [WL[linds[xind, yind], ind5] for xind in xinds, yind in yinds][:]
ws = [W[linds[xind, yind], ind5] for xind in xinds, yind in yinds][:]
writedlm("./out/figures/decay_figure/adapted_basis_function_lvl_5_q_$(q).csv", hcat(xs, ys, vs))
writedlm("./out/figures/decay_figure/basis_function_lvl_5_q_$(q).csv", hcat(xs, ys, ws))

# Plotting a basis function on lvl 6
ind6 = 3600
linds = LinearIndices((n, n))
xinds = (Int(n / 2)) : n
yinds = (Int(n / 2)) : n
xs = [x[1, linds[xind, yind]] for xind in xinds, yind in yinds][:]
ys = [x[2, linds[xind, yind]] for xind in xinds, yind in yinds][:]
vs = [WL[linds[xind, yind], ind6] for xind in xinds, yind in yinds][:]
ws = [W[linds[xind, yind], ind6] for xind in xinds, yind in yinds][:]
writedlm("./out/figures/decay_figure/adapted_basis_function_lvl_6_q_$(q).csv", hcat(xs, ys, vs))
writedlm("./out/figures/decay_figure/basis_function_lvl_6_q_$(q).csv", hcat(xs, ys, ws))

surface(xs, ys, vs)







# scene = Scene()
# fig = Figure()
# #heatmap(log10.(1e-7 .+ abs.(L'[:, end : -1 : 1001])), colormap=["white", "silver"])
# # heatmap(log10.(1e-7 .+ abs.(L'[:, end : -1 : 1001])), colormap=["white", "silver"])
# heatmap(log10.(1e-7 .+ abs.(L'[:, end : -1 : 1])), colormap=["white", rust])
# tightlimits!(fig.axis)

lind = LinearIndices((2 ^ 4, 2 ^ 4))
r = 4
c = 4^4 .+ (vec([lind[i, j] for i = 1 : r: 2^4, j = 1 : r : 2^4]) .- 1) * 3 .+ rand(1:3, length(1 : r : 2^4)^2) .+ 63

xinds = 1 : n
yinds = 1 : n
xs = [x[1, linds[xind, yind]] for xind in xinds, yind in yinds][:]
ys = [x[2, linds[xind, yind]] for xind in xinds, yind in yinds][:]
vs = [sum(WL[linds[xind, yind], c]) for xind in xinds, yind in yinds][:]
ws = [sum(W[linds[xind, yind], c]) for xind in xinds, yind in yinds][:]

writedlm("./out/figures/decay_figure/multicolor_adapted_basis_function_lvl_5_q_$(q).csv", hcat(xs, ys, vs))
writedlm("./out/figures/decay_figure/multicolor_basis_function_lvl_5_q_$(q).csv", hcat(xs, ys, ws))
surface(xs, ys, vs)


xcs = Float64[]
ycs = Float64[]
vcs = Float64[]

xregs = Float64[]
yregs = Float64[]
vregs = Float64[]

n_rows = 4^5
# n_rows = N
cart_inds = CartesianIndices((n_rows, length(c)))
for (ind, val) in enumerate(Lnorm[1 : n_rows, c])
    xind = cart_inds[ind][2] 
    yind = cart_inds[ind][1] 
    push!(xcs, xind)
    push!(ycs, yind)
    push!(vcs, val)
end



fig = Figure()
fig[1, 1] = ax = Axis(fig)
ax.yreversed=true
ax.aspect=AxisAspect(1)


# heatmap!(ax, xcs, ycs, log10.(1e-10 .+ abs.(vcs)), colormap=["white", silver])
heatmap!(ax, xcs, ycs, abs.(vcs) .≥ 1 * 1e-3, colormap=["white", joshua], interpolate=false)

ylims!(ax, (minimum(c), maximum(c)))
ax.yreversed=true

vlines!(ax, 0.5 : 1 : 16.5, color=silver)
# scatter!(ax, xregs, yregs, colormap=["white", silver])
# scatter!(ax, xcs, ycs, colormap=["white", joshua])

writedlm("./out/figures/decay_figure/color_log10L_q_$(q).csv", hcat(xcs, ycs, vcs))
writedlm("./out/figures/decay_figure/regular_log10L_q_$(q).csv", hcat(xregs, yregs, vregs))

save("./out/figures/decay_figure/color_heatmap_log10L.pdf", fig, resolution = (1000, 1000))


# Plotting a basis function on lvl 3
ind_bl = 19
linds = LinearIndices((n, n))
xinds = 1 : Int(n / 2)
yinds = 1 : Int(n / 2)
# xinds = 1 : n
# yinds = 1 : n
xs = [x[1, linds[xind, yind]] for xind in xinds, yind in yinds][:]
ys = [x[2, linds[xind, yind]] for xind in xinds, yind in yinds][:]
vs = [WL[linds[xind, yind], ind_bl] for xind in xinds, yind in yinds][:]
ws = [W[linds[xind, yind], ind_bl] for xind in xinds, yind in yinds][:]


writedlm("./out/figures/decay_figure/scatter_bl.csv", hcat(xs, ys, vs))

ind_tl = 28
linds = LinearIndices((n, n))
xinds = 1 : Int(n / 2)
yinds = (Int(n / 2)) : n
# xinds = 1 : n
# yinds = 1 : n
xs = [x[1, linds[xind, yind]] for xind in xinds, yind in yinds][:]
ys = [x[2, linds[xind, yind]] for xind in xinds, yind in yinds][:]
vs = [WL[linds[xind, yind], ind_tl] for xind in xinds, yind in yinds][:]
ws = [W[linds[xind, yind], ind_tl] for xind in xinds, yind in yinds][:]



writedlm("./out/figures/decay_figure/scatter_tl.csv", hcat(xs, ys, vs))

ind_br = 55
linds = LinearIndices((n, n))
xinds = (Int(n / 2) ) : n
yinds = 1 : Int(n / 2)
# xinds = 1 : n
# yinds = 1 : n
xs = [x[1, linds[xind, yind]] for xind in xinds, yind in yinds][:]
ys = [x[2, linds[xind, yind]] for xind in xinds, yind in yinds][:]
vs = [WL[linds[xind, yind], ind_br] for xind in xinds, yind in yinds][:]
ws = [W[linds[xind, yind], ind_br] for xind in xinds, yind in yinds][:]


writedlm("./out/figures/decay_figure/scatter_br.csv", hcat(xs, ys, vs))

ind_tr = 64
linds = LinearIndices((n, n))
xinds = (Int(n / 2)) : n
yinds = (Int(n / 2)) : n
# xinds = 1 : n
# yinds = 1 : n
xs = [x[1, linds[xind, yind]] for xind in xinds, yind in yinds][:]
ys = [x[2, linds[xind, yind]] for xind in xinds, yind in yinds][:]
vs = [WL[linds[xind, yind], ind_tr] for xind in xinds, yind in yinds][:]
ws = [W[linds[xind, yind], ind_tr] for xind in xinds, yind in yinds][:]

writedlm("./out/figures/decay_figure/scatter_tr.csv", hcat(xs, ys, vs))

# writedlm("./out/figures/decay_figure/adapted_basis_function_lvl_6_q_$(q).csv", hcat(xs, ys, vs))
