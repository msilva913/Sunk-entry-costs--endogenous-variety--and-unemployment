## Simulation  and calculation of moments
#=
Here we follow standard practice and Coles and Kelishomi 2018 by
1) Generating monthly series
2) Converting to quarterly
3) Applying HP filter (lam=100,000)
-> can consider other filters/growth rates, HP filter induces spurious autocorrelations
=#
using Serialization
using PrettyTables
include("time_series_fun.jl")

model_output = deserialize("model_output_GS.jls")
model, targets, par, sol = model_output 

flag_IR = false
flag_logdev = false #express results in LEVELS
T_SM = 100_000
moments_vars = [:u, :v, :θ, :z]
# sim_SM = simulate_model(model, sol_mat, T_SM, eta, SS, flag_IR, flag_logdev)


# # Multiply by 100 to 
# sim_data = 100 .*DataFrame(sim_SM, varnames)
# # Extract variable symbols to be used for calculating moments
# 

# sim_data = sim_data[!, moments_vars]

# # Plot in levels 
# p = Plots.plot(layout=(2, 3), size=(1000,600), 
# legend=true, alpha=0.6)

# # Top left plot
# plot!(p[1], sim_data.u, label=L"u", subplot=1)
# plot!(p[2], sim_data.v, label=L"v")
# plot!(p[3], sim_data.θ, label=L"θ")
# plot!(p[4], sim_data.z, label=L"z")
# #Plots.savefig("simulated_data_levels.pdf")
# display(p)

#######################################################
# Express results in log deviations
flag_logdev = true 
sim_SM = simulate_model(model, sol_mat, T_SM, eta, SS, flag_IR, flag_logdev)
sim_data = DataFrame(sim_SM, varnames)
sim_data = sim_data[!, moments_vars]
# Levels @. exp(sim_data.u)*ss.u


#moments(sim_data, :z, [:z]; lags =2, verbose=true)

# Convert monthly data to quarterly 
sim_data_q = monthly_to_quarterly(sim_data)

# HP and Hamilton filters

sim_data_hp = copy(sim_data_q)
#sim_data_ham = copy(sim_data_q)
#sim_data_growth = copy(sim_data_q)
for x in moments_vars
    sim_data_hp[!, x] .= hp_filter(sim_data_q[!, x], 100_000)
    #sim_data_ham[!, x] .= hamilton_filter(sim_data_q[!, x])
    #sim_data_growth[!, x] .= growth_filter(sim_data_q[!, x])
end

# Calculate moments
    # Set up correlations as Shimer 2005 (w/o job finding rate): 
@show mom = moments(sim_data_hp, :z, [:z]; lags=2)
mom = mom[:,1:4]
pretty_table(mom, backend = Val(:latex))