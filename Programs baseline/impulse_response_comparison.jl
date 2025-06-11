

include("impulse_response_plots.jl")

irf_df = deserialize("irf_z.jls")
irf_df_gen_CES = deserialize("irf_z_gen_CES.jls")


# 1) Compare z shock: baseline and no variety effects
gen_irf_comp(irf_z_gen_CES, irf_df, ["No variety effects", "Baseline"])
savefig("variety_effects_irf_z_comparison.pdf")


# 2) Compare δ shock to CK
irf_δ = deserialize("irf_δ.jls")
irf_δ_CK = deserialize("irf_δ_CK.jls")
gen_irf_comp_simp(irf_δ, irf_δ_CK, ["Baseline", "Coles and Kelishomi"])
savefig("CK_comparison.png")

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
