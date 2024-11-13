# (c) Alvaro Salazar-Perez and Hernán D. Seoane
# "Perturbating and Estimating DSGE models in Julia
# This version 2023

# If you use these codes, please cite us.
# These codes are free and without any guarantee. Use them at your own risk.
# If you find typos, please let us know.

#v1.7+
using MKL
using PyPlot
using DataFrames
#v1.7- 
#BLAS.vendor() 
#:mkl

include("solution_functions.jl")
include("steady_state.jl")


## Model
# Adjustments
    flag_order      = 1
    flag_deviation  = true
    flag_SSsolver   = false

# Parameters
    @syms f_e s zbar δbar b ϕ ρ σ ε A η_L F κ ξ_inv ρ_z σ_z ρ_δ σ_δ
    parameters      = [f_e; s; zbar; δbar; b; ϕ; ρ; σ; ε; A; η_L; F; κ; ξ_inv; ρ_z; σ_z; ρ_δ; σ_δ]
    estimate        = []
    position        = []
    priors          = (;)

        # Transformations
        β = 1/(1+ρ)
        ξ = 1/ξ_inv
        τbar = 1 - (1-δbar)*(1-s)
        μ = ε/(ε-1)

# Variables
    @syms  z δ u θ q K L u v v_pret e K Q  N p N_e ν_f d_f w_R w L_e L_c Y_c C λ Y labor_prod
    @syms zp  δp θp qp Kp Lp up vp v_pretp ep Kp Qp Np pp N_ep ν_fp d_fp w_Rp wp L_ep L_cp Y_cp Cp λp Yp labor_prod_p
    
    
    x               = [u; N; v_pret; z; δ] # predetermined
    y               = [θ; q; L; v; e; K; Q; p; N_e; ν_f; d_f; w_R; w; L_e; L_c; Y_c; C; λ; Y; labor_prod]
    xp              = [up; Np; v_pretp; zp; δp]
    yp              = [θp; qp; Lp; vp; ep; Kp; Qp; pp; N_ep; ν_fp; d_fp; w_Rp; wp; L_ep; L_cp; Y_cp; Cp; λp; Yp; labor_prod_p]
    variables       = [x; y; xp; yp]
    varnames = vcat(Symbol.(x), Symbol.(y))

# Shock
    @syms epsilon
    ex               = [epsilon]
    eta = Array([0.0; 0.0; 0.0; -σ_z; σ_δ]) # size of state space

        
    nx = length(x) 
    ny = length(y)
    nvar = nx + ny 
    ne = length(ex) 

    # Array of symbolics storing model equtions 
    f = fill(Sym("x"), nvar)
    # Equilibrium conditions
        # Job creation condition -> θ
        f[1]  =  κ + K/q - β*λp/λ*(1-δ)*(w_Rp-wp-Kp+(1-s)*(κ+Kp/qp))
        # Marginal revenue product -> w_R
        f[2] = w_R - p*z*zbar/μ
        # Wage equation -> w
        f[3] = w - (ϕ*(w_R-K+θ/(1-δ)*(K+q*κ)) +(1-ϕ)*b)
        # Value of a vacancy -> Q
        f[4] = Q - (e/F)^(ξ_inv)
        # Expected discounted difference in vacancy value -> K
        f[5] = K - (Q-β*λp/λ*(1-δ)*Qp)
        # Market tightness -> v
        f[6] = θ - v/u 
        # Vacancy filling probability -> q
        f[7] = q - A*θ^(-η_L) 
        # Aggregate labor -> L
        f[8] = L - (1-u)
        # Composition of labor -> L_c
        f[9] = L - (L_c+L_e) 
        # Relative price -> p
        f[10] = p - N^(1/(ε-1))
        # Lagrangian multiplier -> λ
        f[11] = λ - C^(-σ)
        # Firm value -> ν_f
        f[12] = ν_f - β*(1-δ)*λp/λ*(ν_fp+d_fp)
        # Retail output: resources -> Y_c
        f[13] = Y_c - p*z*zbar*L_c
        # Retail output: expenditure -> C
        f[14] = Y_c - (C+F/(1+ξ_inv)*(e/F)^(1+ξ_inv)+κ*v*q)
        # Business entrants -> N_e
        f[15] = N_e - L_e*z*zbar/f_e 
        # Firm value relative to price 
        f[16] = ν_f - p*f_e/μ 
        # Output = expenditure
        f[17] = Y - (Y_c+ν_f*N_e)
        # Output = income 
        f[18] = Y - (w_R*L+N*d_f)
        # LOM of vacancies 
        f[19] = v - (v_pret + e)
        # Predetermined vacancies 
        f[20] = v_pretp - (1-δ)*((1-q)*v+s*(1-u))
        # LOM of unemployment
        f[21] = up - ((1-(1-δ)*(θ*q))*u + τ*(1-u))
        # LOM of firms 
        f[22] = Np - (1-δ)*(N+N_e)
        # Profits 
        #f[23] = d_f - Y_c/(N*ε)
        # Data-consistent labor productivity
        f[23] = labor_prod - Y/(p*L)

        # Exogenous processes
        f[24]  =   log(zp) - ρ_z * log(z)
        f[25] =    log(δp) -  ρ_δ * log(δ)

    # Steady state     
        # Values 

        # Initial parameters: targets and normalizations
        N_s = 1.0
        w_s = 1.0
        z_s = 1.0
        fbar = 0.41
        qbar = 0.8
        x_v = 0.20 
        labor_share = 0.66

        p_s = N_s^(1/(ε-1))
        N_es = δ/(1-δ)*N_s
        q_s = qbar/(1-δ)
        θ_s = fbar/qbar
        u_s = τ/(τ+(1-δ)*(θ_s*q_s))
        L_s = 1 - u_s 
        v_s = θ_s*u_s
        e_s = δ*(v_s+1-u_s)
        v_prets = v_s - e_s 
        recruiter_share= (δ+(ρ+δ)*(ε-1))/(δ+(ρ+δ)*ε)
        w_wR = labor_share/recruiter_share
        w_Rs = w_s/(w_wR)
        surplus_ratio = (ρ+τ)/(1-δ)*(1/(q_s*x_v))
        K_s = (w_Rs-w_s)/(1+surplus_ratio) 
        κ   = (1-x_v)/x_v*K_s/q_s
        f_e = (μ-1)*zbar*(L_s/N_s)*(1-δ)/(δ*μ+ρ)
        ν_fs =p_s*f_e/μ
        d_fs = (ρ+δ)/(1-δ)*ν_fs
        L_es = N_es*f_e/zbar
        L_cs = L_s - L_es 
        Y_cs = p_s*zbar*L_cs
        Q_s = K_s*(1+ρ)/(ρ+δ)
        C_s = Y_cs -  F/(1+ξ_inv)*(e_s/F)^(1+ξ_inv) -  κ*v_s*q_s 
        λ_s = C_s^(-σ)
        Y_s = Y_cs + ν_fs*N_es
        labor_prod_s = Y_s/(p_s*L_s)
        #x               = [u; N; v_pret; z] # predetermined
        #y               = [θ; q; L; v; e; K; Q; p; N_e; ν_f; d_f; w_R; w; L_e; L_c; Y_c; C; λ; Y; labor_prod]
        # Vector
        SS_block  = [log(x) for x in [u_s, N_s, v_prets, z_s, δ_s, θ_s, q_s, L_s, v_s, e_s, K_s, Q_s, p_s, N_es, ν_fs, d_fs, w_Rs, w_s, L_es, L_cs, Y_cs, C_s, λ_s, Y_s, labor_prod_s]]
        # vertical concatenate: represent both current and future variables
        SS = vcat(SS_block, SS_block)

        #PAR_SS = [ALPHA; BETA; DELTA; RHO; SIGMA; MUU; AA]
        PAR_SS = parameters[:]


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
    targets = (labor_share=labor_share, dest_ann=0.06, r_ann=0.04, f =fbar, η_L=0.6, q=qbar, sep=0.031, b_ratio=0.71, 
            x_v=0.20, ξ_inv=1, ε=4.3, σ=1.5, N=N_s, w=w_s)
    cal = calibrate_labor_share(targets)
    
    
    # Shock values
    ρ_z = 0.975
    σ_z = 0.007

@unpack  f_e, δ, s, z, b, ϕ, ρ, σ, ε, A, η_L, F, κ, ξ_inv = cal

zbar = z 
δbar = δ
PAR     =   [f_e; s; zbar; δbar; b; ϕ; ρ; σ; ε; A; η_L; F; κ; ξ_inv; ρ_z; σ_z; ρ_δ; σ_δ ]

sol = solution_interface(model, PAR)
@unpack ss, SS, sol_mat, eta = sol

## Simulation  and calculation of moments
#=
Here we follow standard practice and Coles and Kelishomi 2018 by
1) Generating monthly series
2) Converting to quarterly
3) Applying HP filter (lam=100,000)
-> can consider other filters/growth rates, HP filter induces spurious autocorrelations
=#




flag_IR = false
flag_logdev = true #express results in log deviations
T_SM = 100_000
sim_SM = simulate_model(model, sol_mat, T_SM, eta, SS, flag_IR, flag_logdev)
# Multiply by 100 to express results in percentage deviations
sim_data = 100 .*DataFrame(sim_SM, varnames)
# Extract variable symbols to be used for calculating moments
moments_vars = [:u, :v, :θ, :z, :δ ]

sim_data = sim_data[!, moments_vars]
#moments(sim_data, :z, [:z]; lags =2, verbose=true)

# Convert monthly data to quarterly 
sim_data_q = monthly_to_quarterly(sim_data)

# HP and Hamilton filters

function filter_data(data, type=:hp, λ=100_000)
sim_data_hp = copy(sim_data)
sim_data_ham = copy(sim_data)
sim_data_growth = copy(sim_data)
for x in moments_vars
    sim_data_hp[!, x] .= hp_filter(sim_data[!, x], 1600)
    sim_data_ham[!, x] .= hamilton_filter(sim_data[!, x])
    sim_data_growth[!, x] .= growth_filter(sim_data[!, x])
end

# Calculate moments
    # Set up correlations as Shimer 2005 (w/o job finding rate): 
@show mom = moments(sim_data_growth, :z, [:v, :θ, :δ, :z]; lags=2)
#moments(sim_data_ham, :z, [:z, :v]; lags =2, verbose=true)

# Calculate moments 

#sim_dat
#using Plots
#Plots.plot(sim_SM)



