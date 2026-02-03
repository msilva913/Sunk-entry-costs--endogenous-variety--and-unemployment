include("steady_state.jl")

############################################################
# --- Calibration Targets and Steady-State Calculation -----
############################################################

# Define calibration targets (edit these as needed)
targets = (
    X_Y=0.015,        # recruiting cost share of output
    Xc_Y=0.20,         # fixed cost share of output (Abraham, Bormans, Konings, Roeger)
    dest_ann=0.0754,   # annual product destruction rate
    dest_end_frac=0.5, # endogenous share of destruction rate
    f=0.41,            # job-finding rate, 
    η_L=0.6,           # elasticity of matching fun wrt unemployment
    q=0.8,             # vacancy filling rate,
    sep=0.031,         # aggregate separation rate , 
    b_ratio=0.71,      # ratio of unemployment benefits to wage,
    x_v=1.0, 
    ξ_inv=1, 
    r_ann=0.04,        # annual discount rate
    ε=4.3,             # Elasticity of substitution (BGM, Compustat)
    σ=1.0,             # Inverse IES
    N=1.0,             # SS mass of forms (normalization)
    w=1.0,             # SS wage (normalization)
    #ψ=1.5)
)


# Calibrate and compute steady state
cal = calibrate_shares(targets)
steady = steady_state(cal)

############################################################
# ----------- Unpack Calibrated Variables ------------------
############################################################

@unpack θ, δ_e, x_c, ρ, L_c, L_e, w, w_int, L, N, N_e, K, Q, q, f, ν_f, d_f, C, Y_c, Y, X_v, X, X_c,
         labor_share, sunk_vac_cost_share, vacancy_share, x_v, M, e, v, u, π_s,
         profit_share_rec, profit_share_ret = steady

@unpack f_e, δ, z, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, x_m, s, f_m, ψ = cal

# Composites
ψ_c = ψ/(ψ+1)
τ = 1-(1-δ_e)*(1-s)
μ = ε/(ε-1)


############################################################
# -------------------- Consistency Checks ------------------
############################################################

# Separation elasticity ψ F(x_c)/(1-F(x_c))
dest_elast_1 = dest_elast(cal, x_c)
@show surv_prob = (x_c/f_m)^ψ
@show dest_elast_2 = ψ*surv_prob/(1-surv_prob)
@assert abs(dest_elast_1 - dest_elast_2) < 1e-12

# Normalizations and target matches
@assert abs(steady.N - targets.N) < 1e-12
@assert abs(steady.f*(1-δ_e) - targets.f) < 1e-12
@assert abs(steady.ρ - 1.0) < 1e-12
@assert abs(targets.X_Y - vacancy_share) < 1e-12
@assert abs(targets.Xc_Y - X_c/Y) < 1e-12


labor_share_alt = 1.0 - profit_share_rec - profit_share_ret
@assert abs(labor_share - labor_share_alt) < 1e-12

# Retail profits 
@assert abs(N*d_f - π_s*Y_c) < 1e-12
 - Y_c/ε - X_c
@assert abs(N*d_f - Y_c/ε + X_c) < 1e-12

# Wage bill share check
#@assert abs(w_int*L/Y - (1+X/Y)*(δ+(ρ+δ)*(ε-1))/(δ+(ρ+δ)*ε)) < 1e-12

# Output identity
@assert abs(Y-C-ν_f*N_e) < 1e-12

# Production and cost identities
@assert (Y_c - ρ*z*L_c) < 1e-12
@assert (Y_c - C - X - X_c) < 1e-12
#@assert abs(Y_c+ν_f*N_e - ρ*z*L_c - ρ*z*L_e/μ) < 1e-12

@assert abs(ρ*z*L_e/μ-ν_f*N_e) < 1e-12

# Vacancy cost share
@assert abs(κ/(κ+K/q) -(1-x_v)) < 1e-12

# Wage equation/JCC
@assert abs((w_int-w-K - (r+τ)/(1-δ_e)*(κ+K/q))) < 1e-12
@assert abs(w - ϕ*(w_int-K+ θ*(K+q* κ))- (1-ϕ)*b) < 1e-12

# Entrant consistency
@assert abs(v - ((1-δ_e)*((1-q)*v+s*(1-u))+e)) < 1e-12
@assert abs(e - δ_e*(v+1-u)) < 1e-12
@assert abs(N_e - L_e*z/f_e) < 1e-12

# Surplus ratio and capital check
surplus_ratio = (r+τ)/(1-δ_e)*(1/(q*x_v))
@assert abs(K - (1-ϕ)/ϕ*(w-b)/(surplus_ratio +  θ*(1/x_v))) < 1e-12

# Consistency x_m with Q
@assert abs(Q - e^ξ_inv*x_m) < 1e-12

# Cutoff
x_c - ρ*f_e/μ*((r+δ_e)/((1-δ_e)*(μ*π_s))*(μ-1) + 1)
@assert abs(x_c - ρ*f_e/μ*((r+δ_e)/(1-δ_e)*(r+δ_e+ψ_c*(1-δ_e))/((1-ψ_c)*(r+δ_e))+1))<1e-12

############################################################
# ----------- Display Key Steady-State Ratios --------------
############################################################

# TO DO: Curves to be revised
@show N_jcc(steady.θ, cal)
@show N_res(steady.θ, cal)
@show labor_share
@show C/Y
@show X/Y
@show sunk_vac_cost_share
@show vacancy_share
@show X/(q*v*w)
@show x_v

# Profits
@show ((w_int-w)*L-X)/Y
@show (N*d_f)/Y

# Value of unemployment benefit
@show b/w_int
@show b/(Y/L)
@show b/z

@show M/(12*Y)

############################################################
# ----------- Further Checks and Experiments ---------------
############################################################

# Check consistency of calibrate function (with ϕ added)
"""
targets_aug = (targets..., ϕ=cal.ϕ)
cal_aug = calibrate(targets_aug)
steady2 = steady_state(cal_aug)

# Study implication of low b
targets_lowb = (targets..., b_ratio=0.41)
cal_lowb = calibrate_labor_share(targets_lowb)
steady_state(cal_lowb)
# For very low b, labor share can be matched; for high b, may need negative bargaining power.

# Entry elasticity ξ → ∞ (ξ_inv → 0)
para = (cal..., ξ_inv=0.00)
steady = steady_state(para)
@assert abs(steady.Q-1.0) < 1e-12
@assert abs(steady.X_v-steady.e) < 1e-12
@assert abs(steady.K - (cal.ρ+cal.δ)/(1+cal.ρ)) < 1e-12

# Alternate calibration: direct ϕ specification (example, commented)
# targets_ϕ = (ϕ=0.6, dest_ann=0.06, r_ann=0.04, f=0.41, η_L=0.6, q=0.8, sep=0.031, b_ratio=0.71, ξ_inv=1, ε=4, σ=1.5, N=1, w=1.0)
# cal_ϕ = calibrate(targets_ϕ)
# steady_ϕ = steady_state(cal_ϕ)

# Very high elasticity of substitution
para_highε = (para..., ε=1e12)
steady_highε = steady_state(para_highε)
@assert abs(steady_highε.Q-1.0) < 1e-12
@assert abs(steady_highε.X_v-steady_highε.e) < 1e-12
@assert abs(steady_highε.K - (cal.ρ+cal.δ)/(1+cal.ρ)) < 1e-12

# Risk aversion check
targets_risk = (
    labor_share = 0.66, dest_ann = 0.10, r_ann = 0.04, f = 0.41, η_L = 0.6, q = 0.8,
    sep = 0.031, b_ratio = 0.71, x_v = 0.1, ξ_inv = 1/0.265, ε = 4.3, σ = 0.0, N = 1.0, w = 1.0
)
cal_risk = calibrate_labor_share(targets_risk)
steady_risk = steady_state(cal_risk)
@show steady_risk.labor_share
@show steady_risk.sunk_vac_cost_share
@show steady_risk.x_v
@show steady_risk.C/steady_risk.Y
@show steady_risk.M/(12*steady_risk.Y)
"""