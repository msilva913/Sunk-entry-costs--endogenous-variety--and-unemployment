include("steady_state_CK.jl")

targets = (ϕ=0.6, r_ann=0.04, f =0.41, η_L=0.6, q=0.8, sep=0.0304, b_ratio=0.71, ξ_inv=1/0.265, w=1.0)
cal = calibrate(targets)

steady = steady_state(cal)
@unpack θ, w, L, K, q, C, Y, X, labor_share, M = steady
@unpack δ, z, b, ϕ, ρ, A, η_L, ξ_inv, F = cal
# Accuracy checks
@assert abs(steady.f*(1-cal.δ) - targets.f) < 1e-12



@assert abs(Y - C - X) < 1e-12
@assert abs(Y - z*L) < 1e-12
@assert abs(w - 1.0) < 1e-12



# Steady-state ratios
@show labor_share
@show C/Y
@show M/(12*Y)

cal_table = calibration_table(cal, targets)



################################################################
