# Based on original code by Alvaro Salazar-Perez and Hernán D. Seoane
# Modified by Mario Silva

using DataFrames
using Serialization
cd(@__DIR__)
#v1.7- 
#BLAS.vendor() 
#:mkl

include("solution_functions.jl")
include("steady_state_CK.jl")
include("impulse_response_plots_CK.jl")
include("time_series_fun.jl")

function solution_interface(model, PAR)
    eta     =   eval_ShockVAR(PAR)
    PAR_SS  =   eval_PAR_SS(PAR)
    SS      =   eval_SS(PAR_SS)
    SS_err  =   eval_SS_error(PAR_SS, SS)
    deriv   =   eval_deriv(PAR_SS, SS)
    SS_max = maximum(abs.(SS_err))
    argmax(abs.(SS_err))
    println("Residuals: $SS_max")

    ss = NamedTuple(zip(model.varnames, exp.(SS[1:model.nvar])))
    sol_mat = solve_model(model, deriv, eta)
    println("Model solved")
    out = (SS=SS, ss=ss, eta=eta, deriv=deriv, sol_mat=sol_mat)
    return out
end

## Model
# Adjustments
    flag_order      = 1
    flag_deviation  = true
    flag_SSsolver   = false

# Parameters
    @syms  zbar δbar b ϕ ρ A η_L F x_m ξ_inv ρ_z σ_z ρ_δ σ_δ
    parameters      = [zbar; δbar; b; ϕ; ρ; A; η_L; F; x_m; ξ_inv; ρ_z; σ_z; ρ_δ; σ_δ ]
    estimate        = []
    position        = []
    priors          = (;)

    # Transformations
    β = 1/(1+ρ)
    ξ = 1/ξ_inv

# Variables
@syms  z δ u θ q K L v v_pret e K Q  w C Y
@syms  zp δp up θp qp Kp Lp vp v_pretp ep Kp Qp  wp Cp Yp


x               = [u; v_pret; z; δ] # predetermined
y               = [θ; q; L; v; e; K; Q;  w; C; Y]
xp              = [up; v_pretp; zp; δp]
yp              = [θp; qp; Lp; vp; ep; Kp; Qp; wp; Cp; Yp]
variables       = [x; y; xp; yp]
varnames = vcat(Symbol.(x), Symbol.(y))

# Shock
@syms epsilon
ex               = [epsilon]
eta = Array([0.0; 0.0; -σ_z; σ_δ]) # size of state space

    
nx = length(x) 
ny = length(y)
nvar = nx + ny 
ne = length(ex) 

# Array of symbolics storing model equtions 
f = fill(Sym("x"), nvar)
# Equilibrium conditions
    # Job creation condition -> θ
    f[1]  =   K/q - β*(1-δbar*δ)*(zp-wp-Kp +Kp/qp)
    # Wage equation -> w
    f[2] = w - (ϕ*(z-K+θ*K) +(1-ϕ)*b)
    # Value of a vacancy -> Q
    f[3] = Q - (e/F)^(ξ_inv)*x_m
    # Expected discounted difference in vacancy value -> K
    f[4] = K - (Q-β*(1-δbar*δ)*Qp)
    # Market tightness -> v
    f[5] = θ - v/u 
    # Vacancy filling probability -> q
    f[6] = q - A*θ^(-η_L) 
    # Aggregate labor -> L
    f[7] = L - (1-u)
    # Output: Technology
    f[8] = Y - z*L
    # Output -> C
    f[9] = Y - (C+e/(1+ξ_inv)*Q)
    # Business entrants -> N_e
    # LOM of vacancies 
    f[10] = v - (v_pret + e)
    # Predetermined vacancies  -> v_pret
    f[11] = v_pretp - (1-δbar*δ)*((1-q)*v)
    # LOM of unemployment
    f[12] = up - ((1-(1-δbar*δ)*(θ*q))*u +δ*δbar*(1-u))

    # Exogenous processes
    f[13]  =   log(zp) - ρ_z * log(z)
    f[14] =    log(δp) -  ρ_δ * log(δ)

# Steady state   
"""
If known, provide the symbolic expresion of the steady state in 
terms of the parameters.
It must be provided as a column vector, in which the variables 
should be ordered in the following way: [x; y; xp; yp].
If the steady state is unknown leave it as an empty vector (SS=[]), 
so that the program tries to estimate it.
"""  
function SS_symbolics(parameters::Vector{Sym{PyObject}}, targets)

     zbar, δbar, b, ϕ, ρ, A, η_L, F, x_m, ξ_inv, ρ_z, σ_z, ρ_δ, σ_δ = parameters
    # Initial parameters: targets and normalizations/ leave parameters as symbolic to be populated with calibration
    @unpack  ϕ, f, q, w = targets
    w_s = w 
    fbar = f 
    qbar = q 

    q_s = qbar/(1-δbar)
    θ_s = fbar/qbar
    u_s =  δbar/(δbar+(1-δbar)*(θ_s*q_s))
    L_s = 1 - u_s 
    v_s = θ_s*u_s
    e_s = δbar*(v_s+1-u_s)
    v_prets = v_s - e_s 

    # Rescale zbar to be consistent with wage=1
    surplus_ratio = (ρ+δbar)/(1-δbar)*(1/(q_s))
    K_s = (1-ϕ)/ϕ*(w_s-b)/(surplus_ratio +  θ_s)
    zbar = surplus_ratio*K_s + w_s + K_s

    Y_s = zbar*L_s
    Q_s = K_s*(1+ρ)/(ρ+δbar)
    C_s = Y_s - e_s/(1+ξ_inv)*Q_s

    # Shocks
    z_s = 1.0 
    δ_s = 1.0
    s_s = 1.0

    #x               = [u; N; v_pret; z] # predetermined
    #y               = [θ; q; L; v; e; K; Q; p; N_e; ν_f; d_f; w_int; w; L_e; L_c; Y_c; C; λ; Y; labor_prod]
    # Vector
    SS_block  = [log(x) for x in [u_s, v_prets, z_s, δ_s, θ_s, q_s, L_s, v_s, e_s, K_s, Q_s, w_s, C_s, Y_s]]
    # vertical concatenate: represent both current and future variables
    SS = vcat(SS_block, SS_block)
    return SS
end

#PAR_SS = [ALPHA; BETA; DELTA; RHO; SIGMA; MUU; AA]
PAR_SS = parameters[:]

targets = (ϕ=0.6, r_ann=0.04, f =0.41, η_L=0.6, q=0.8, sep=0.034, b_ratio=0.7, ξ_inv=1/0.265, w=1.0)
SS = SS_symbolics(parameters, targets)

# Values 
cal = calibrate(targets)



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
#parameters     = [zbar; δbar; b, ϕ, ρ, A, η_L, F, ξ_inv, ρ_z, σ_z, ρ_δ, σ_δ]


# Shock values (from Coles and Kelishomi), monthly frequency
ρ_z = 0.965
σ_z = 0.007
ρ_δ = 0.875
σ_δ = 0.042

@unpack  z, δ, b, ϕ, ρ, A, η_L, F, ξ_inv = cal

zbar = z 
δbar = δ
PAR     =   [zbar; δbar; b; ϕ; ρ; A; η_L; F; ξ_inv; ρ_z; σ_z; ρ_δ; σ_δ]

sol = solution_interface(model, PAR)
@unpack ss, SS, sol_mat, eta = sol
pprint(ss)

## Checks 
ss2 = steady_state(cal)


# Export: model, targets, PAR, sol (save output using serialization)
model_output_CK = (model, targets, PAR, sol)
#serialize("model_output.jls", model_output)


## Impulse responses ##
eta_z = zero(eta) # Tech shock
eta_δ = zero(eta) # Product destruction shock

eta_z[3] = eta[3]
eta_δ[4] = eta[4]

flag_IR = true
flag_logdev = true
T_IR = 120 # 10 years

# Destruction rate shock: consistent with Beveridge curve
irf_δ= simulate_model(model, sol_mat, T_IR, eta_δ, SS, flag_IR, flag_logdev) 
irf_δ = 100 .*DataFrame(irf_δ, varnames)

irf_δ_red = irf_δ[1:60, :]
gen_irf(irf_δ_red)
savefig("δ_shock_CK.png")
#Plots.savefig("dest_shock.pdf")
serialize("irf_δ_CK.jls", irf_δ)

"""
# Technology shock
irf_z= simulate_model(model, sol_mat, T_IR, eta_z, SS, flag_IR, flag_logdev)
irf_z = 100 .*DataFrame(irf_z, varnames)
gen_irf(irf_z)
Plots.savefig("z_shock.pdf")
savefig("z_shock.png")
#serialize("irf_z.jls", irf_z)
"""

