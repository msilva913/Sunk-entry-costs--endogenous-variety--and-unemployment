
include("run_solution_core.jl")
# Values 

#PAR_SS = [ALPHA; BETA; DELTA; RHO; SIGMA; MUU; AA]
PAR_SS = parameters[:]

targets = (targets..., x_v=1.0)

SS = SS_symbolics(parameters, targets)
cal = calibrate_labor_share(targets)


f = gen_model_equations()
# Procesing the model (no adjustment needed)                
model = (parameters = parameters, estimate = estimate, estimation = position,
        npar = length(parameters), ns = length(estimate), 
        priors = priors,
        x = x, y = y, xp = xp, yp = yp, variables = variables,
        varnames=varnames, #store symbols of variable names
        nx = nx, ny = ny, nvar = nvar,
        e = ex, eta = eta,
        ne = ne,
        f = f,
        nf = nvar,
        SS = SS, PAR_SS = PAR_SS,
        flag_order = flag_order, flag_deviation = flag_deviation, flag_SSsolver = flag_SSsolver)
process_model(model)


## Solution
# Parametrization: need to handle dependent parameters 
#parameters      = [f_e; δ; s; zbar; b; ϕ; ρ; σ; ε; A; η_L; F; κ; ξ_inv; ρ_z; μ_z]

# Shock values (from Coles and Kelishomi), monthly frequency
ρ_z = 0.965
σ_z = 0.007
ρ_δ = 0.875
σ_δ = 0.042
ρ_s = 0.875
σ_s = 0.0042

@unpack  f_e, δ, s, z, b, ϕ, ρ, σ, ε, A, η_L, F, κ, ξ_inv = cal

zbar = z 
δbar = δ
sbar = s
PAR     =   [f_e; zbar; δbar; sbar; b; ϕ; ρ; σ; ε; A; η_L; F; κ; ξ_inv; ρ_z; σ_z; ρ_δ; σ_δ; ρ_s; σ_s ]

sol = solution_interface(model, PAR)
@unpack ss, SS, sol_mat, eta = sol
# Export: model, targets, PAR, sol (save output using serialization)
model_output = (model, targets, PAR, sol)
serialize("model_output.jls", model_output)


## Impulse responses ##
eta_z = zero(eta) # Tech shock
eta_δ = zero(eta) # Product destruction shock
eta_s = zero(eta) # Idiosyncratic separation shock
eta_τ = zero(eta) # Common separation shock

eta_z[4] = eta[4]
eta_δ[5] = eta[5]
eta_s[6] = eta[6]

eta_τ[5] = eta[5]
eta_τ[6] = eta[6]
flag_IR = true
flag_logdev = true
T_IR = 60 # 10 years


# Technology shock
irf_z= simulate_model(model, sol_mat, T_IR, eta_z, SS, flag_IR, flag_logdev)
irf_z = 100 .*DataFrame(irf_z, varnames)
gen_irf(irf_z)
#savefig("z_shock.png")
serialize("irf_z_kappa.jls", irf_z)

# Destruction rate shock: consistent with Beveridge curve
irf_δ= simulate_model(model, sol_mat, T_IR, eta_δ, SS, flag_IR, flag_logdev) 
irf_δ = 100 .*DataFrame(irf_δ, varnames)
serialize("irf_δ_kappa.jls", irf_δ)
gen_irf(irf_δ)
#Plots.savefig("dest_shock.pdf")

"""
# Idiosyncratic job separation shock
irf_s= simulate_model(model, sol_mat, T_IR, eta_s, SS, flag_IR, flag_logdev) 
irf_s = 100 .*DataFrame(irf_s, varnames)
gen_irf(irf_s)
Plots.savefig("s_shock.pdf")

#Commmon separation shock
irf_τ= simulate_model(model, sol_mat, T_IR, eta_τ, SS, flag_IR, flag_logdev) 
irf_τ = 100 .*DataFrame(irf_τ, varnames)
gen_irf(irf_τ)
Plots.savefig("common_shock.pdf")
"""


