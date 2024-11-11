

# Comparison to high epsi calibration 
targets2 = (labor_share=labor_share, dest_ann=0.06, r_ann=0.04, f =fbar, η_L=0.6, q=qbar, sep=0.031, b_ratio=0.71, 
x_v=0.20, ξ_inv=1, ε=100, σ=1.5, N=N_s, w=w_s)
cal2 = calibrate_labor_share(targets2)

@unpack  f_e, δ, s, z, b, ϕ, ρ, σ, ε, A, η_L, F, κ, ξ_inv = cal2

zbar = z 
δbar = δ
PAR2     =   [f_e; s; zbar; δbar; b; ϕ; ρ; σ; ε; A; η_L; F; κ; ξ_inv; ρ_z; σ_z; ρ_δ; σ_δ ]
sol2 = solution_interface(model, PAR2)
sol_mat2 = sol2.sol_mat
SS2 = sol2.SS

sim_IR2 = simulate_model(model, sol_mat2, T_IR, eta_z, SS2, flag_IR, flag_logdev)
irf_df2 = 100 .*DataFrame(sim_IR2, varnames)

gen_irf_comp(irf_df, irf_df2, ["Baseline", " ε=100"])

# Comparison to high δ calibration: we maintain aggregate worker separations at 3.1%
targets3 = (labor_share=labor_share, dest_ann=0.2, r_ann=0.04, f =fbar, η_L=0.6, q=qbar, sep=0.031, b_ratio=0.71, 
x_v=0.20, ξ_inv=1, ε=4.3, σ=1.5, N=N_s, w=w_s)
cal3 = calibrate_labor_share(targets3)

@unpack  f_e, δ, s, z, b, ϕ, ρ, σ, ε, A, η_L, F, κ, ξ_inv = cal3

zbar = z 
δbar = δ
PAR3     =   [f_e; s; zbar; δbar; b; ϕ; ρ; σ; ε; A; η_L; F; κ; ξ_inv; ρ_z; σ_z; ρ_δ; σ_δ ]
sol3 = solution_interface(model, PAR3)
sol_mat3 = sol3.sol_mat
SS3 = sol3.SS
sim_IR3 = simulate_model(model, sol_mat3, T_IR, eta_z, SS3, flag_IR, flag_logdev)
irf_df3 = 100 .*DataFrame(sim_IR3, varnames)

gen_irf_comp(irf_df, irf_df3, ["Baseline", "High δ"])
# Compare irf's for different δ
# Compare irf's for s vs. δ
# Compare irf's for different ξ_inv
# How to decompose aggregate separations--which can be computed using Shimer's approach on unemployment flows--into series for s and δ
# Compute moments for filtered data (HP and Hamilton)
# function to compute moments 
# filter