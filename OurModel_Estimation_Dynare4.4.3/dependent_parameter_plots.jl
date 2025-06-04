
using Parameters, CSV, StatsBase, Statistics, Random
using NLsolve
using DataFrames
using PrettyPrinting
using LaTeXStrings
cd("C:/Users/msilva913/Documents/GitHub/Sunk_entry_costs_endogenous_variety_unemployment/Programs baseline/")
include(pwd()*"/steady_state.jl")
include(pwd()*"/generate_tables.jl")

# Trace plots of parameters 
cd(@__DIR__)
Traceplot_parameter = matopen("Traceplot_parameter.mat")
Traceplot_parameter = read(Traceplot_parameter, "Traceplot_parameter")
# For steady-state-based dependent parameters, we can omit Shocks
Traceplot_parameter = Traceplot_parameter[:, 1:(end-6)]
col_names = [:σ, :b, :x_v, :ξ_inv, :δ, :ε]
df = DataFrame(Traceplot_parameter, :auto)
rename!(df, col_names)

#labor_share, dest_ann, r_ann, f, η_L, q, sep, b_ratio, x_v, ξ_inv, ε, σ, N, w = targets
#targets = (labor_share=0.66, dest_ann=0.10, r_ann=0.04, f =0.41, η_L=0.6, q=0.8, sep=0.031, b_ratio=0.71, x_v=0.1, ξ_inv=1/0.265, ε=4.3, σ=1.0, N=1.0, w=1.0)
n = size(df, 1)
targets_df = DataFrame(
labor_share = 0.66*ones(n),
dest_ann = 1.0 .-((1. .-df[!,:δ]).^12),
r_ann = 0.04*ones(n),
f = 0.41*ones(n),
η_L = 0.6*ones(n),
q = 0.8*ones(n),
sep = 0.031*ones(n), 
b_ratio = df[!, :b],
x_v = df[!, :x_v],
ξ_inv = df[!, :ξ_inv],
ε = df[!, :ε],
σ = df[!, :σ],
N = 1.0*ones(n),
w = 1.0*ones(n),
)

results = [calibrate_labor_share(targets_df[row,:]) for row in 1:n]

dep_parameters_df = DataFrame(results)

κ = dep_parameters_df[!, :κ]
ϕ = dep_parameters_df[!,:ϕ]
using KernelDensity
filtered_data(data) = filter(!isnan, data)

density = kde(filtered_data(ϕ))

fig, ax = plt.subplots()
ax.plot(density.x, density.density, color="orange", alpha=0.3)
ax.fill_between(density.x, density.density, color="gold", alpha=0.3, zorder=1)
ax.set_xlabel("ϕ", fontsize=12)
ax.set_ylabel("Density", fontsize=12)
display(fig)


posterior_mean= matopen("posterior_mean.mat")
posterior_mean = read(posterior_mean, "posterior_mean")

delta_mean = posterior_mean["delta"]
delta_ann_mean = 1 - (1-delta_mean)^12

targets = (labor_share=0.66, 
           dest_ann=delta_ann_mean, 
           r_ann=0.04, 
           f =0.41, 
           η_L=0.6, #elast. of matching function
           q=0.8, 
           sep=0.031, 
           b_ratio=posterior_mean["b_ratio"], 
            x_v=posterior_mean["xi_inv"], 
            ξ_inv=posterior_mean["xi_inv"], # congestion elasticity
            ε=posterior_mean["epsi"], # elasticity of sub.
            σ=posterior_mean["sigma"], # log utility
            N=1, w=1.0)

cal = calibrate_labor_share(targets)
ss = steady_state(cal)
cal_table = calibration_table(cal, targets)

output = IOBuffer()
show(output, MIME("text/latex"),cal_table)
cal_table_tex = String(take!(output))
print(cal_table_tex)

df_shares = shares_table(ss)

output = IOBuffer()
show(output, MIME("text/latex"),df_shares)
df_shares_table = String(take!(output))
print(df_shares_table)
