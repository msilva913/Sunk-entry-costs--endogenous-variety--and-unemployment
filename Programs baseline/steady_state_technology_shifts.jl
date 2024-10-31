using PyPlot
using Parameters, CSV, StatsBase, Statistics, Random
using NLsolve
using Roots, Optim, LeastSquaresOptim
using PrettyPrinting
using Distributions
cd(@__DIR__)
include("steady_state.jl")

function θ_prod(para)
    zgrid = range(para.z, stop=para.z*0.9, length=10)
    θ_grid = similar(zgrid)
    u_grid = similar(zgrid)
    N_grid = similar(zgrid)
    labor_prod = similar(zgrid)
    for (i,z) in enumerate(zgrid)
        para = (para..., z=z)
        steady = steady_state(para)
        θ_grid[i] = steady.θ
        u_grid[i] = steady.u 
        N_grid[i] = steady.N 
        labor_prod[i] = steady.labor_prod
    end
    elast = (θ_grid[2]-θ_grid[1])/(zgrid[2]-zgrid[1])*zgrid[1]/θ_grid[1]
    return zgrid, labor_prod, θ_grid, u_grid, N_grid, elast
end

targets = (labor_share=0.66, dest_ann=0.06, r_ann=0.04, f =0.41, η_L=0.6, q=0.8, sep=0.031, b_ratio=0.71, x_v=0.2, ξ_inv=1, ε=4, σ=1.5, N=1, w=1)
cal = calibrate_labor_share(targets)
steady = steady_state(cal)
zgrid, labor_prod, θ_grid, u_grid, N_grid, elast = θ_prod(cal)

targets2 = (targets..., ε=100)
cal2 = calibrate_labor_share(targets2)
zgrid2, labor_prod2, θ_grid2, u_grid2, N_grid2, elast2 = θ_prod(cal2)

targets3 = (targets..., x_v=1.0)
#targets3 = (targets..., σ=0.0)
cal3 = calibrate_labor_share(targets3)
zgrid3, labor_prod3, θ_grid3, u_grid3, N_grid3 = θ_prod(cal3)


fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
θ_vars = [θ_grid, θ_grid2, θ_grid3]
u_vars = [u_grid, u_grid2, u_grid3]
N_vars = [N_grid, N_grid2, N_grid3]
vars = [θ_vars, u_vars, N_vars]
labels = ["Tightness", "Unemployment", "Number of firms"]
for i in 1:3
    var = vars[i]
    ax[i].plot(1:10, var[1], label="ε=4, x_v=0.2")
    ax[i].plot(1:10, var[2], label="ε=100, x_v=0.2")
    ax[i].plot(1:10, var[3], label="ε=4, x_v=1.0")
    #ax[i].plot(1:10, var[3], label="ε=4, σ=0.0")
    ax[i].set_xlabel("Periods")
    ax[i].set_title(labels[i])
    ax[i].legend()
end
plt.savefig("response_tightness_markup.pdf")
display(fig)

# fig, ax = plt.subplots()
# ax.plot(1:10, labor_prod, label="ε=4")
# ax.plot(1:10, labor_prod2, label="ε=10")
# ax.set_xlabel("Technology")
# ax.legend()
# display(fig)
