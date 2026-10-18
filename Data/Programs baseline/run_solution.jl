
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
# Shock processes: SHOCKS_BASE in run_solution_core.jl (part6b; replaced by the Q2 VAR).
# Never assign numbers to the symbolic ρ_z, σ_z, … (gen_model_equations reads them).
SHK = SHOCKS_BASE
PAR = build_PAR(cal, SHK)

sol = solution_interface(model, PAR, SS_numeric(PAR, targets))   # never the symbolic SS (M10)
@unpack ss, SS, sol_mat = sol   # do NOT unpack eta: sol.eta is Matrix{Sym} (symbolic zeros), not Float64
# Export: model, targets, PAR, sol (save output using serialization)
model_output = (model, targets, PAR, sol)
serialize("model_output.jls", model_output)


## Impulse responses ##
# Single-shock columns are built by state name (shock_col), since R8 changed the state vector.
flag_IR     = true
flag_logdev = true
T_IR        = 60   # 5 years at monthly frequency

eta_z_col = shock_col(:z, SHK.σz)
eta_δ_col = shock_col(:δ, SHK.σδ)
eta_s_col = shock_col(:s, SHK.σs)

# ne=1: one shock fired at a time; sim_shocks has 1 column, matching eta_*_col
model_1shock = (; model..., ne = 1)

# Technology shock
irf_z = simulate_model(model_1shock, sol_mat, T_IR, eta_z_col, SS, flag_IR, flag_logdev)
irf_z = 100 .* DataFrame(irf_z, varnames)
gen_irf(irf_z)
savefig("z_shock.png")
serialize("irf_z.jls", irf_z)

# Permanent exit (δ) shock: core Beveridge curve prediction
irf_δ = simulate_model(model_1shock, sol_mat, T_IR, eta_δ_col, SS, flag_IR, flag_logdev)
irf_δ = 100 .* DataFrame(irf_δ, varnames)
gen_irf(irf_δ)
savefig("delta_shock.png")
serialize("irf_delta.jls", irf_δ)

# Worker separation (s) shock
irf_s = simulate_model(model_1shock, sol_mat, T_IR, eta_s_col, SS, flag_IR, flag_logdev)
irf_s = 100 .* DataFrame(irf_s, varnames)
gen_irf(irf_s)
savefig("s_shock.png")
serialize("irf_s.jls", irf_s)


