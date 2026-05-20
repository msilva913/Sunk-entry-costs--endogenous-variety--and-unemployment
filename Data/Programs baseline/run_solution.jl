
include("run_solution_core.jl")

PAR_SS = parameters[:]
SS     = SS_symbolics(parameters, targets)
cal    = calibrate_shares(targets)

f = gen_model_equations()
model = (parameters = parameters, estimate = estimate, estimation = position,
         npar      = length(parameters), ns = length(estimate),
         priors    = priors,
         x = x, y = y, xp = xp, yp = yp, variables = variables,
         varnames  = varnames,
         nx = nx, ny = ny, nvar = nvar,
         e = ex, eta = eta,
         ne = ne,
         f = f,
         nf = nvar,
         SS = SS, PAR_SS = PAR_SS,
         flag_order = flag_order, flag_deviation = flag_deviation,
         flag_SSsolver = flag_SSsolver)
process_model(model)


## Solution
# ── Shock process parameters (monthly frequency) ──────────────────────────────
# Source: part6b_shock_persistence_cyclical.py
#   Bivariate VAR(1) on HP-filtered (λ=1600) log(z) and log(δ),
#   Cholesky-identified with z ordered first (δ ⊥ z by construction).
#   Quarterly parameters converted to monthly via ρ_m = ρ_q^(1/3) and
#   unconditional-variance-matching for σ.  Sample: 1992Q3–2019Q4 (N=110).
#   δ shock excludes the endogenous z→δ channel; that pathway lives in
#   model equilibrium equations (exit threshold x_c).
#   s shock: no direct data counterpart after separating δ and z — treated
#   as a residual; set to small value pending SMM estimation.

ρ_z = 0.902    # quarterly: 0.735; monthly via ρ_q^(1/3)
σ_z = 0.0092   # structural z innovation SD (log-cycle units, monthly)
ρ_δ = 0.592    # quarterly: 0.207; monthly via ρ_q^(1/3)
σ_δ = 0.0669   # structural δ innovation SD ⊥ z (log-cycle units, monthly)
ρ_s = 0.90     # placeholder — to be estimated via SMM
σ_s = 0.010    # placeholder — to be estimated via SMM

# Values at posterior mode (uncomment when estimation is wired up):
# ρ_z = posterior_mode["rho_z"]
# σ_z = posterior_mode["sigma_z"]
# ρ_δ = posterior_mode["rho_delta"]
# σ_δ = posterior_mode["sigma_delta"]
# ρ_s = posterior_mode["rho_s"]
# σ_s = posterior_mode["sigma_s"]

# ── Parameter vector ──────────────────────────────────────────────────────────
# Order must match @syms declaration in run_solution_core.jl
@unpack f_e, δ, s, z, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, x_m, ψ, f_m, p_0 = cal

zbar = z
δbar = δ
sbar = s

PAR = [f_e; zbar; δbar; sbar; b; ϕ; r; σ; ε; A; η_L; κ; ξ_inv; x_m; ψ; f_m; p_0;
       ρ_z; σ_z; ρ_δ; σ_δ; ρ_s; σ_s]

sol = solution_interface(model, PAR)
@unpack ss, SS, sol_mat, eta = sol
# Export: model, targets, PAR, sol (save output using serialization)
model_output = (model, targets, PAR, sol)
serialize("model_output.jls", model_output)


## Impulse responses ##
# eta is 6×3: rows = state variables [u,N,v_pret,z,δ,s], cols = [z,δ,s] shocks
flag_IR     = true
flag_logdev = true
T_IR        = 60   # 5 years at monthly frequency

# Each IRF fires one column of eta at t=1
eta_z = eta[:, 1]   # z shock (col 1)
eta_δ = eta[:, 2]   # δ shock (col 2)
eta_s = eta[:, 3]   # s shock (col 3)

# Technology shock
irf_z = simulate_model(model, sol_mat, T_IR, eta_z, SS, flag_IR, flag_logdev)
irf_z = 100 .* DataFrame(irf_z, varnames)
gen_irf(irf_z)
savefig("z_shock.png")
serialize("irf_z.jls", irf_z)

# Permanent exit (δ) shock: core Beveridge curve prediction
irf_δ = simulate_model(model, sol_mat, T_IR, eta_δ, SS, flag_IR, flag_logdev)
irf_δ = 100 .* DataFrame(irf_δ, varnames)
gen_irf(irf_δ)
savefig("delta_shock.png")
serialize("irf_delta.jls", irf_δ)

# Worker separation (s) shock
irf_s = simulate_model(model, sol_mat, T_IR, eta_s, SS, flag_IR, flag_logdev)
irf_s = 100 .* DataFrame(irf_s, varnames)
gen_irf(irf_s)
savefig("s_shock.png")
serialize("irf_s.jls", irf_s)


