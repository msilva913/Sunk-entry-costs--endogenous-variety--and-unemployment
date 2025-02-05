using Plots
import Plots:plot, savefig

cd(@__DIR__)
include("steady_state_general_CES.jl")
Plots.gr()
default(linewidth=2, grid=true, fontfamily="Computer Modern")

targets = (labor_share=0.66, dest_ann=0.06, r_ann=0.04, f =0.41, η_L=0.6, q=0.8, sep=0.031, b_ratio=0.71, 
    x_v=0.2, ξ_inv=1, ε=4.3, σ=1.0, N=1, w=1,  ζ=0)
cal = calibrate_labor_share(targets)
steady = steady_state(cal)


T = 100
θ_grid = range(0.1, stop=1.0, length=T)
N_res_grid = N_res.(θ_grid, Ref(cal))
#N_jcc_grid = N_jcc.(θ_grid, Ref(cal))
θ_jcc(cal)
   
# Initial plot of $\theta$ and $N$
# Initial plot
#plt = plot(θ_grid, N_jcc_grid, label="N: job creation condition", linewidth=2)  # Adjust the alpha value here
plt=plot(θ_grid, N_res_grid, label="N: aggregate resource constraint", linewidth=2)  # Adjust the alpha value here
xlabel!(L"θ")
ylabel!(L"N")
vline!([θ_jcc(cal)], line=:dash, color=:black, label="Job creation condition "* L"θ", linewidth=2)  # Adjust the alpha value here
hline!([steady.N], line=:dash, color=:black, label="Steady "*"N", linewidth=0.5)  # Adjust the alpha value here
plot!(legend=true)
savefig("equilibrium_curves_gen_CES.pdf")
display(plt)

# Check for second Equilibrium
#steady = steady_state(cal, init=0.01)

#N_res is increasing with markup (decreasing with ε)-> incentive to place more resources in business formation
para = (cal..., ε=cal.ε*0.9)
N_res_grid2 = N_res.(θ_grid, Ref(para))
p = plot( size=(600, 400), legend=:right)
plot!(p[1], θ_grid, N_res_grid, label="N: resource constraint curve: ε=$(cal.ε)", line=:solid, color="blue")
plot!(p[1], θ_grid, N_res_grid2, label="N: resource constraint curve: ε=$(para.ε)", line=:dash, color="red")
vline!([θ_jcc(cal)], line=:solid, label="Job creation condition: ε=$(cal.ε)", linewidth=2, color="blue")  # Adjust the alpha value here
vline!([θ_jcc(para)], line=:dash, label="Job creation condition: ε=$(para.ε)", linewidth=2, color="red")  # Adjust the alpha value here
xlabel!(p[1], L"θ")
ylabel!(p[1], L"N")
savefig("curves_epsi_shift_general_CES.pdf")
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

