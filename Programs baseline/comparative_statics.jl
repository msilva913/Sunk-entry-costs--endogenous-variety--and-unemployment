using Plots
import Plots:plot, savefig

cd(@__DIR__)
include("steady_state.jl")
Plots.gr()
default(linewidth=2, grid=true, fontfamily="Computer Modern")
targets = (labor_share=0.66, dest_ann=0.10, f =0.41, η_L=0.6, q=0.8, sep=0.031, b_ratio=0.71, 
            x_v=1.0, ξ_inv=1, X_Y=0.015, C_Y=0.80, σ=1.0, N=1.0, w=1.0)
cal = calibrate_shares(targets)
steady = steady_state(cal)

#########################################
T = 100
ξ_inv_series = range(0.95, 1.2, length=T)
X_Y_series = similar(ξ_inv_series)
for (i, ξ_inv) in enumerate(ξ_inv_series)
    para = (;cal..., ξ_inv=ξ_inv) #merge syntax to update para 
    ss = steady_state(para)
    X_Y_series[i] = ss.vacancy_share
end

################
T = 1000
θ_grid = range(1e-100, stop=0.9, length=T)

N_jcc_grid = N_jcc.(θ_grid, Ref(cal))
N_res_grid = N_res.(θ_grid, Ref(cal))
   
# Initial plot of $\theta$ and $N$
# Initial plot
plt = plot(θ_grid, N_jcc_grid, label="Job creation curve", linewidth=2)  # Adjust the alpha value here
plot!(θ_grid, N_res_grid, label="Resource constraint curve", linewidth=2)  # Adjust the alpha value here
xlabel!(L"θ")
ylabel!(L"N")
vline!([steady.θ], line=:dash, color=:black, linewidth=0.5, label=nothing)  # Adjust the alpha value here
hline!([steady.N], line=:dash, color=:black, linewidth=0.5, label=nothing)  # Adjust the alpha value here
plot!(legend=true)
savefig("equilibrium_curves.pdf")
display(plt)

# Check for second Equilibrium
#steady = steady_state(cal, init=0.01)

#N_res is increasing with markup (decreasing with ε)-> incentive to place more resources in business formation
para = (cal..., ε=cal.ε*0.9)
N_res_grid2 = N_res.(θ_grid, Ref(para))
N_jcc_grid2 = N_jcc.(θ_grid, Ref(para))
p = plot(layout=(1, 2), size=(800, 400), legend=:right)
plot!(p[1], θ_grid, N_res_grid, label="N: resource constraint curve: ε=$(cal.ε)")
plot!(p[1], θ_grid, N_res_grid2, label="N: resource constraint curve: ε=$(para.ε)")
xlabel!(p[1], L"θ")
ylabel!(p[1], L"N")

plot!(p[2], θ_grid, N_jcc_grid, label="N: job creation curve: ε=$(cal.ε)")
plot!(p[2], θ_grid, N_jcc_grid2, label="N: job creation curve: ε=$(para.ε)")
xlabel!(p[2], L"θ")
ylabel!(p[2], L"N")
savefig("curves_epsi_shift.pdf")
display(p)

#N_res decreasing with δ (reduces labor resources, increases effective discounting, decreases ratio of N to N_e)

para = (cal..., δ=cal.δ*1.1)
N_res_grid2 = N_res.(θ_grid, Ref(para))

p = plot(θ_grid, N_res_grid, label="N: aggregate resource constraint")
plot!(θ_grid, N_res_grid2, label="N: aggregate resource constraint, higher destruction rate δ")
xlabel!(L"θ")
ylabel!(L"N")
plot!(legend=true)
display(p)


# Baseline exercise: increase in f_e

para = (cal..., f_e=cal.f_e*1.2)
steady2 = steady_state(para)
N_jcc_grid2 = N_jcc.(θ_grid, Ref(para))
N_res_grid2 = N_res.(θ_grid, Ref(para))

# New plot
plt = plot()
plot!(plt, θ_grid, N_jcc_grid, label="N: job creation condition")
plot!(plt, θ_grid, N_res_grid, label="N: aggregate resource constraint")
plot!(plt, θ_grid, N_res_grid2, label="N: aggregate resource constraint, higher sunk entry cost f_e")
vline!(plt, [steady.θ], linestyle=:dash, linecolor=:black, linewidth=0.5, label=false)
hline!(plt, [steady.N], linestyle=:dash, linecolor=:black, linewidth=0.5, label=false)
vline!(plt, [steady2.θ], linestyle=:dash, linecolor=:red, linewidth=0.5, label=false)
hline!(plt, [steady2.N], linestyle=:dash, linecolor=:red, linewidth=0.5, label=false)
plot!(plt, legend=true)
display(plt)
# Increase in F, amount of firms which can create vacancies

para = (cal..., F=cal.F*1.2)
steady2 = steady_state(para)
N_jcc_grid2 = N_jcc.(θ_grid, Ref(para))
N_res_grid2 = N_res.(θ_grid, Ref(para)) # only shifts N_jcc

# New plot
plt = plot()
plot!(plt, θ_grid, N_jcc_grid, label="N: job creation condition")
plot!(plt, θ_grid, N_jcc_grid2, label="N: job creation condition, higher F")
plot!(plt, θ_grid, N_res_grid, label="N: aggregate resource constraint")
vline!(plt, [steady.θ], linestyle=:dash, linecolor=:black, linewidth=0.5, label=false)
hline!(plt, [steady.N], linestyle=:dash, linecolor=:black, linewidth=0.5, label=false)
vline!(plt, [steady2.θ], linestyle=:dash, linecolor=:red, linewidth=0.5, label=false)
hline!(plt, [steady2.N], linestyle=:dash, linecolor=:red, linewidth=0.5, label=false)
plot!(plt, legend=true)



############################
# Numerical differentiation 
function numerical_derivative(f::Function, x::Float64, h::Float64=1e-5)
    return (f(x+h)-f(x-h))/(2*h)
end

function elasticity(f::Function, x::Float64, h::Float64=1e-7)
    f_x = f(x) 
    df_dx = numerical_derivative(f,x,h)
    return (x/f_x)*df_dx
end

function w_fun(z, cal)
   para = (cal..., z=z)
   steady = steady_state(para)
   return steady.w/steady.p
end

function labor_prod_fun(z, cal)
   para = (cal..., z=z)
   steady = steady_state(para)
   return steady.Y/(steady.p*steady.L)
end

targets = (labor_share=0.66, dest_ann=0.10, r_ann=0.04, f =0.41, η_L=0.6, q=0.8, sep=0.031, b_ratio=0.71,
 x_v=1.0, ξ_inv=1.0, ε=4.3, σ=1.0, N=1, w=1)
cal = calibrate_labor_share(targets)
#steady = steady_state(para)
elast_w_z = elasticity(z -> w_fun(z, cal),cal.z)
elast_lp_z = elasticity(labor_prod_fun, cal.z)

