include("steady_state_CK_parallel_calibration.jl")
# Coles and Kelishomi generalized with separate s and delta shocks
############################################################
# --- Calibration Targets and Steady-State Calculation -----
############################################################
targets = (
    X_Y      = 0.015,
    dest_ann = 0.0754,
    f        = 0.41,
    η_L      = 0.6,
    q        = 0.8,
    sep      = 0.031,
    b_ratio  = 0.71,
    ξ_inv    = 1.0,
    r_ann    = 0.04,
    w        = 1.0
)

cal = calibrate_shares(targets)

steady = steady_state(cal)
@unpack θ, w, L, K, q, C, Y, X, entrant_share, vacancy_cost_share, e, u, v, M, profit_share, labor_share = steady
@unpack δ, s, z, b, ϕ, ρ, A, η_L, ξ_inv, x_m, F = cal
# Accuracy checks
@assert abs(steady.f*(1-cal.δ) - targets.f) < 1e-12



@assert abs(Y - C - X) < 1e-12
@assert abs(Y - z*L) < 1e-12
#@assert abs(w - 1.0) < 1e-12



# Steady-state ratios
@show u 
@show v
@show vacancy_cost_share
@show entrant_share 
@show profit_share 
@show M/(12*Y)

cal_table = calibration_table(cal, targets)



################################################################
