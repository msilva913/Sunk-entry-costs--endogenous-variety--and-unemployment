include("steady_state.jl")

############################################################
# --- Calibration Targets and Steady-State Calculation -----
############################################################

# Define calibration targets (edit these as needed)
targets = (
    X_Y      = 0.015,
    dest_ann = 0.0754,
    f        = 0.41,
    η_L      = 0.6,
    q        = 0.8,
    sep      = 0.031,
    b_ratio  = 0.71,
    x_v      = 0.5,
    ξ_inv    = 1.0,
    r_ann    = 0.04,
    ε        = 4.3,
    σ        = 1.0,
    N        = 1.0,
    w        = 1.0
)

# Calibrate and compute steady state
cal = calibrate_shares(targets)
steady = steady_state(cal)

############################################################
# ----------- Unpack Calibrated Variables ------------------
############################################################

@unpack θ, p, L_c, L_e, w, w_int, L, N, N_e, K, Q, q, f, ν_f, d_f, C, Y_c, Y, X_v, X,
         labor_share, sunk_vac_cost_share, vacancy_share, x_v, M, e, v, u = steady

@unpack f_e, τ, δ, z, b, ϕ, ρ, σ, ε, A, η_L, κ, ξ_inv, F, x_m, s = cal

μ = ε/(ε-1)

############################################################
# -------------------- Consistency Checks ------------------
############################################################

# Normalizations and target matches
@assert abs(steady.N - targets.N) < 1e-12
@assert abs(steady.f*(1-cal.δ) - targets.f) < 1e-12
@assert abs(steady.p - 1.0) < 1e-12
@assert abs(targets.X_Y - vacancy_share) < 1e-12

# Labor share check (theoretical formula)
@assert abs(labor_share - (1+X/Y)*(w/w_int)*(δ+(ρ+δ)*(ε-1))/(δ+(ρ+δ)*ε)) < 1e-12

# Profit share of consumption output
@assert abs(N*d_f - Y_c/ε) < 1e-12

# Wage bill share check
@assert abs(w_int*L/Y - (1+X/Y)*(δ+(ρ+δ)*(ε-1))/(δ+(ρ+δ)*ε)) < 1e-12

# Output identity
@assert abs(Y-C-ν_f*N_e) < 1e-12

# Production and cost identities
@assert abs(Y+X - p*z*L_c - p*z*L_e/μ) < 1e-12
@assert abs(p*z*L_c - w_int*L - N*ν_f*ρ/(1-δ)) < 1e-12
@assert abs(p*z*L_e/μ-ν_f*N_e) < 1e-12

# Vacancy cost share
@assert abs(κ/(κ+K/q) -(1-x_v)) < 1e-12

# Wage equation/JCC
@assert abs((w_int-w-K - (ρ+τ)/(1-δ)*(κ+K/q))) < 1e-12
@assert abs(w - ϕ*(w_int-K+ θ*(K+q* κ))- (1-ϕ)*b) < 1e-12

# Entrant consistency
@assert abs(v - ((1-δ)*((1-q)*v+s*(1-u))+e)) < 1e-12
@assert abs(e - δ*(v+1-u)) < 1e-12
@assert abs(N_e - L_e*z/f_e) < 1e-12

# Surplus ratio and capital check
surplus_ratio = (ρ+τ)/(1-δ)*(1/(q*x_v))
@assert abs(K - (1-ϕ)/ϕ*(w-b)/(surplus_ratio +  θ*(1/x_v))) < 1e-12

# Consistency x_m with Q
@assert abs(Q - (e/F)^ξ_inv*x_m) < 1e-12

############################################################
# ----------- Display Key Steady-State Ratios --------------
############################################################

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