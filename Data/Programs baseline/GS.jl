# Based on original code by Alvaro Salazar-Perez and Hernán D. Seoane
# Modified by Mario Silva

using MKL
using DataFrames, Parameters
using Roots
using Serialization
using LaTeXStrings
using PyCall
cd(@__DIR__)
#v1.7- 
#BLAS.vendor() 
#:mkl

function calibrate_GS(targets)
    @unpack zbar, dest_ann, r_ann, f, q, wage_elast, τ, b, ξ_inv, γ = targets

    ρ = (1+r_ann)^(1/12)-1
    β = 1/(1+ρ)
    δ = 1-(1-dest_ann)^(1/12)
    
    # Correct job finding and vacancy filling probablities
    f = f/(1-δ)
    q = q/(1-δ)

    θ = f/q
    u = τ/(τ+(1-δ)*f)
    v = θ*u 

    # Level parameter of matching function 
    A = f/θ^(1-γ)
    s = (τ-δ)/(1-δ)

    e = δ*(v+1-u)

    function loss(ϕ)
        K = q*(1-δ)*(1-ϕ)*(zbar-b)/(ρ+τ+ϕ*f+q*(1-δ)*(1-ϕ))
        w = (ϕ*(zbar-K+θ/(1-δ)*K) +(1-ϕ)*b)
        out = (ϕ*zbar/w - wage_elast)/wage_elast
        return out, K
    end
    
    ϕ = find_zero(x -> loss(x)[1], (0.01, 0.99))
    out, K = loss(ϕ)
    Q = K*(1+ρ)/(ρ+δ) 
    F = e/Q^(1/ξ_inv)

    #zbar δ s b ϕ ρ A γ F ξ_inv
    cal = (zbar=zbar, δ=δ, s=s, b=b, ϕ=ϕ, ρ=ρ,  A=A,  γ=γ, F=F, ξ_inv=ξ_inv)
    return cal
end

include("solution_functions.jl")
#include("steady_state.jl")
include("impulse_response_plots_GS.jl")
#include("time_series_fun.jl")

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
    # dev = zeros(model.nvar)
    # common_keys = intersect(keys(ss), keys(steady))
    # for (i, field) in enumerate(common_keys)
    #     dev[i] = ss[field] -steady[field]
    # end
    # print(maximum(abs.(dev)))
    # @btime sol_mat = solve_model(model, deriv, eta)
    eta = map(float, eta)
    sol_mat = solve_model(model, deriv, eta)
    println("Model solved")
    out = (SS=SS, ss=ss, eta=eta, deriv=deriv, sol_mat=sol_mat)
    return out
end

## Model
# Adjustments
    flag_order      = 2
    flag_deviation  = true
    flag_SSsolver   = false

# Parameters
    @syms zbar δ s b ϕ ρ A γ F ξ_inv ρ_z σ_z
    parameters      = [zbar; δ; s; b; ϕ; ρ; A; γ; F; ξ_inv; ρ_z; σ_z]
    estimate        = []
    position        = []
    priors          = (;)

    # Transformations
    β = 1/(1+ρ)
    ξ = 1/ξ_inv
    τbar = 1 - (1-δ)*(1-s)

# Variables
@syms  u v_pret z θ q v e K Q w
@syms up v_pretp zp θp qp vp ep Kp Qp wp


x               = [u; v_pret; z] # predetermined
y               = [θ; q; v; e; K; Q; w]
xp              = [up; v_pretp; zp]
yp              = [θp; qp; vp; ep; Kp; Qp; wp]
variables       = [x; y; xp; yp]
varnames = vcat(Symbol.(x), Symbol.(y))

# Shock
@syms epsilon
ex               = [epsilon]
eta = Array([0.0; 0.0; -σ_z]) # size of state space

    
nx = length(x) 
ny = length(y)
nvar = nx + ny 
ne = length(ex) 

# Array of symbolics storing model equtions 
f = fill(Sym("x"), nvar)
# Equilibrium conditions
    # Job creation condition -> θ
    f[1]  =  K/q - β*(1-δ)*(zp-wp-Kp+(1-s)*Kp/qp)
    # Wage equation -> w
    f[2] = w - (ϕ*(z-K+θ/(1-δ)*K) +(1-ϕ)*b)
    # Value of a vacancy -> Q
    f[3] = Q - (e/F)^(ξ_inv)
    # Expected discounted difference in vacancy value -> K
    f[4] = K - (Q-β*(1-δ)*Qp)
    # Market tightness -> v
    f[5] = θ - v/u 
    # Vacancy filling probability -> q
    f[6] = q - A*θ^(-γ) 
    # LOM of vacancies 
    f[7] = v - (v_pret + e)
    # Predetermined vacancies 
    f[8] = v_pretp - (1-δ)*((1-q)*v+s*(1-u))
    # LOM of unemployment
    f[9] = up - ((1-(1-δ)*(θ*q))*u + τbar*(1-u))

    # Exogenous processes
    f[10]  =   log(zp) - ρ_z * log(z)

# Steady state   
"""
If known, provide the symbolic expresion of the steady state in 
terms of the parameters.
It must be provided as a column vector, in which the variables 
should be ordered in the following way: [x; y; xp; yp].
If the steady state is unknown leave it as an empty vector (SS=[]), 
so that the program tries to estimate it.
"""  
# Values 
#zbar δ s b ϕ ρ A γ F ξ_inv ρ_z σ_z

function SS_symbolics(parameters::Vector{Sym{PyObject}}, targets)

    zbar, δ, s, b, ϕ, ρ, A, γ, F, ξ_inv, ρ_z, σ_z = parameters
    # Initial parameters: targets and normalizations/ leave parameters as symbolic to be populated with calibration
    @unpack f, q  = targets
    fbar = f 
    qbar = q 
    τbar = 1 - (1-δ)*(1-s)

    z_s = zbar;
    q_s = qbar/(1-δ)
    f_s = fbar/(1-δ)
    θ_s = fbar/qbar
    u_s = τbar/(τbar+(1-δ)*(θ_s*q_s))
    v_s = θ_s*u_s
    e_s = δ*(v_s+1-u_s)
    v_pret_s = v_s - e_s 
    Q_s = (e_s/F)^(ξ_inv)

    K_s = q_s*(1-δ)*(1-ϕ)*(z_s-b)/(ρ+τbar+ϕ*f_s+q_s*(1-δ)*(1-ϕ))
    w_s = ϕ*(z_s-K_s+θ_s/(1-δ)*K_s) +(1-ϕ)*b
    #x  = [u; v_pret; z] # predetermined
    #y  = [θ; q; v; e; K; Q; w]
    # Vector
    SS_block  = [log(x) for x in [u_s, v_pret_s, z_s, θ_s, q_s, v_s, e_s, K_s, Q_s, w_s]]
    # vertical concatenate: represent both current and future variables
    SS = vcat(SS_block, SS_block)
    return SS
end

#PAR_SS = [ALPHA; BETA; DELTA; RHO; SIGMA; MUU; AA]
PAR_SS = parameters[:]

targets = (zbar=1.0,
           dest_ann=0.06,
            r_ann=0.04, 
            f=1/2.2, 
            q=1-(1-1/3)^4,
            #ϕ=0.566, 
            wage_elast=0.6,
            τ=0.034, 
            b=0.9, 
            ξ_inv=1.0,
            #ξ_inv=1.0/0.265,
            γ=0.6)   

SS = SS_symbolics(parameters, targets)
cal = calibrate_GS(targets)



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
#parameters      = [f_e; δ; s; zbar; b; ϕ; ρ; σ; ε; A; γ; F; κ; ξ_inv; ρ_z; μ_z]

# Shock values (from Coles and Kelishomi), monthly frequency
ρ_z = 0.979
σ_z = 0.007


@unpack zbar, δ, s, b, ϕ, ρ, A, γ,F, ξ_inv = cal

PAR     =   [zbar; δ; s; b; ϕ; ρ; A; γ; F; ξ_inv; ρ_z; σ_z ]

sol = solution_interface(model, PAR)
@unpack ss, SS, sol_mat, eta = sol
# Export: model, targets, PAR, sol (save output using serialization)
model_output = (model, targets, PAR, sol)
serialize("model_output_GS.jls", model_output)


## Impulse responses ##
flag_IR = true
flag_logdev = true
T_IR = 120 # 10 years


# Technology shock
irf_z= simulate_model(model, sol_mat, T_IR, eta, SS, flag_IR, flag_logdev)
irf_z = 100 .*DataFrame(irf_z, varnames)
gen_irf(irf_z)




