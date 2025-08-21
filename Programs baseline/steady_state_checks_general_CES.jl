include("steady_state_general_CES.jl")

targets = (labor_share=0.66, dest_ann=0.10, r_ann=0.04, f =0.41, η_L=0.6, q=0.8, sep=0.031, b_ratio=0.71, 
x_v=0.1, ξ_inv=1/0.265, ε=4.3, σ=1.0, N=1.0, w=1.0, ζ=0.0)
cal = calibrate_labor_share(targets)

steady = steady_state(cal)
@unpack θ, p, L_c, L_e, w, w_int, L, N, N_e, K, q, ν_f, d_f, C, Y, Y_c, X_v, X, labor_share, sunk_vac_cost_share, x_v, M = steady
@unpack f_e, τ, δ, z, b, ϕ, ρ, σ, ε, A, η_L, κ, ξ_inv, F, s = cal
μ = ε/(ε-1)
# Accuracy checks
@assert p ≈ 1.0
@assert abs(steady.N - targets.N) < 1e-12
@assert abs(steady.f*(1-cal.δ) - targets.f) < 1e-12


@assert abs(N*d_f - Y_c/ε) < 1e-12 # profit share of consumption output
@assert abs(w_int*L/Y - (δ+(ρ+δ)*(ε-1))/(δ+(ρ+δ)*ε)) < 1e-12
@assert abs(Y - w_int*L - N*d_f) < 1e-12
@assert abs(Y - p*z*L_c - p*z*L_e/μ) < 1e-12
@assert abs(p*z*L_c - w_int*L - N*ν_f*ρ/(1-δ)) < 1e-12
@assert abs(p*z*L_e/μ-ν_f*N_e) < 1e-12
@assert abs(labor_share - targets.labor_share) < 1e-12
@assert abs(labor_share - (w/w_int)*(δ+(ρ+δ)*(ε-1))/(δ+(ρ+δ)*ε)) < 1e-12
@assert abs(κ/(κ+K/q) -(1-x_v)) < 1e-12
X_v/X

@assert abs(N_e - L_e*z/f_e) < 1e-12
surplus_ratio = (ρ+τ)/(1-δ)*(1/(q*x_v))
@assert abs(K - (1-ϕ)/ϕ*(w-b)/(surplus_ratio +  θ/(1-δ)*(1/x_v))) < 1e-12

# Steady-state ratios
@show labor_share
@show sunk_vac_cost_share
@show x_v
@show C/Y
@show M/(12*Y)

# Consistency checks: should replicate steady state
#N_jcc(steady.θ, cal)
θ_jcc(cal)
N_res(steady.θ, cal)


################################################################
# Check consistency of calibrate function 
targets = (targets...,  ϕ=cal.ϕ)
cal = calibrate(targets)
steady2 = steady_state(cal)
