

using Serialization
cd(@__DIR__)
include("impulse_response_plots.jl")

irf_z = deserialize("irf_z.jls")
irf_z_gen_CES = deserialize("irf_z_gen_CES.jls")
T=60
#####################
# Main impulse response: 
gen_irf(irf_z[1:T,:])
# calculate impulse response
@show mean(cumsum(irf_z.w_R[2:24])./cumsum(irf_z.labor_prod[2:24]))

gen_irf(irf_δ[1:T,:])
@show mean(cumsum(irf_δ.w_R[2:24])./cumsum(irf_δ.labor_prod[2:24]))

# 1) Compare z shock: baseline and no variety effects
gen_irf_comp(irf_z_gen_CES[1:T,:], irf_z[1:T,:], ["No variety effects", "Baseline"])
savefig("variety_effects_irf_z_comparison.png")
savefig("variety_effects_irf_z_comparison.pdf")


# 2) Compare δ shock to CK (assuming $δ shock accounts for all separtions)
irf_δ_alt = deserialize("irf_δ_alt.jls")
irf_δ_CK = deserialize("irf_δ_CK.jls")
gen_irf_comp_simp(irf_δ_alt, irf_δ_CK[1:T,:], ["Baseline", "Coles and Kelishomi"])
savefig("CK_comparison.png")
savefig("CK_comparison.pdf")


# 3) Examine role of risk neutrality 
irf_δ = deserialize("irf_δ.jls")
irf_δ_risk_neutral = deserialize("irf_δ_risk_neutral.jls")
gen_irf_comp(irf_δ, irf_δ_risk_neutral, ["Baseline", "Risk neutral"])
savefig("risk_aversion_comparison_delta_shock.png")
savefig("risk_aversion_comparison_delta_shock.pdf")


irf_z_risk_neutral = deserialize("irf_z_risk_neutral.jls")
gen_irf_comp(irf_z[1:T,:], irf_z_risk_neutral[1:T,:], ["Baseline", "Risk neutral"])
savefig("risk_aversion_comparison_z_shock.png")
savefig("risk_aversion_comparison_z_shock.pdf")


# 4) #Examine role of highly elastic vacancy creation
# Less elastic vacancy creation amplifies separation shocks but dampens tech shocks.
irf_z_elastic = deserialize("irf_z_elastic.jls")
irf_δ_elastic = deserialize("irf_δ_elastic.jls")
gen_irf_comp(irf_δ, irf_δ_elastic, ["Baseline", "ξ_inv=0.5"])
savefig("elastic_comparison_delta_shock.png")
savefig("elastic_comparison_delta_shock.pdf")
gen_irf_comp(irf_z, irf_z_elastic, ["Baseline", "ξ_inv=0.5"])
savefig("elastic_comparison_tech_shock.png")
savefig("elastic_comparison_tech_shock.pdf")


# 5) Compare with high b
irf_z_b = deserialize("irf_z_b.jls")
irf_δ_b = deserialize("irf_δ_b.jls")
gen_irf_comp(irf_z, irf_z_b, ["Baseline", "b=0.92"])
savefig("b_comparison_tech_shock.png")
savefig("b_comparison_tech_shock.pdf")
gen_irf_comp(irf_δ, irf_δ_b, ["Baseline", "b=0.92"])
savefig("b_comparison_delta_shock.png")
savefig("b_comparison_delta_shock.pdf")

# 6) Compare with lower kappa
irf_z_kappa = deserialize("irf_z_kappa.jls")
irf_δ_kappa = deserialize("irf_δ_kappa.jls")
gen_irf_comp(irf_z, irf_z_kappa, ["Baseline", "κ=0"])
savefig("kappa_comparison_tech_shock.png")
savefig("kappa_comparison_tech_shock.pdf")
gen_irf_comp(irf_δ, irf_δ_kappa, ["Baseline", "κ=0"])
savefig("kappa_comparison_delta_shock.png")
savefig("kappa_comparison_delta_shock.pdf")
##################################################################
"""
out = deserialize("model_output.jls")
model, targets, PAR = out
sol = solution_interface(model, PAR)
@unpack ss, SS, sol_mat, eta = sol
# Focus on technology shocks
eta_z = zero(eta) # Tech shock
eta_z[4] = eta[4]

targets2 = (targets..., ε=100.0 )
cal2 = calibrate_labor_share(targets2)
@unpack  f_e, δ, s, z, b, ϕ, ρ, σ, ε, A, η_L, F, κ, ξ_inv = cal2
PAR2     =   [f_e; z; δ; s; b; ϕ; ρ; σ; ε; A; η_L; F; κ; ξ_inv; ρ_z; σ_z; ρ_δ;
 σ_δ; ρ_s; σ_s ]
sol2 = solution_interface(model, PAR2)
sol_mat2 = sol2.sol_mat
SS2 = sol2.SS

flag_IR = true
flag_logdev = true
sim_IR2 = simulate_model(model, sol_mat2, T_IR, eta_z, SS2, flag_IR, flag_logdev)
irf_df2 = 100 .*DataFrame(sim_IR2, varnames)
gen_irf_comp(irf_df, irf_df2, ["Baseline", " ε=100"])
Plots.savefig("irf_comp_epsi.pdf")


# High δ calibration: we maintain aggregate worker separations at 3.1%
targets3 = (targets..., dest_ann=0.2)
cal3 = calibrate_labor_share(targets3)

@unpack  f_e, δ, s, z, b, ϕ, ρ, σ, ε, A, η_L, F, κ, ξ_inv = cal3
PAR3     =   [f_e; z; δ; s; b; ϕ; ρ; σ; ε; A; η_L; F; κ; ξ_inv; ρ_z; σ_z; ρ_δ;
 σ_δ; ρ_s; σ_s ]

sol3 = solution_interface(model, PAR3)
sol_mat3 = sol3.sol_mat
SS3 = sol3.SS
sim_IR3 = simulate_model(model, sol_mat3, T_IR, eta_z, SS3, flag_IR, flag_logdev)
irf_df3 = 100 .*DataFrame(sim_IR3, varnames)

gen_irf_comp(irf_df, irf_df3, ["Baseline", "High δ"])
Plots.savefig("irf_comp_delta.pdf")


# Compare irf's for different δ
# Compare irf's for s vs. δ
# Compare irf's for different ξ_inv
# How to decompose aggregate separations--which can be computed using Shimer's approach on unemployment flows--into series for s and δ
# Compute moments for filtered data (HP and Hamilton)
# function to compute moments 
# filter
"""