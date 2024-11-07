# (c) Alvaro Salazar-Perez and Hernán D. Seoane
# "Perturbating and Estimating DSGE models in Julia
# This version 2023

# If you use these codes, please cite us.
# These codes are free and without any guarantee. Use them at your own risk.
# If you find typos, please let us know.

cd(@__DIR__)
using MKL
using PyPlot
using DataFrames


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
        @syms  z δ u θ q K L u v v_pret e K Q  N p N_e ν_f d_f w_R w L_e L_c Y_c C λ Y 
        @syms zp  δp θp qp Kp Lp up vp v_pretp ep Kp Qp Np pp N_ep ν_fp d_fp w_Rp wp L_ep L_cp Y_cp Cp λp Yp 
       
        
        x               = [u; N; v_pret; z; δ] # predetermined
        y               = [θ; q; L; v; e; K; Q; p; N_e; ν_f; d_f; w_R; w; L_e; L_c; Y_c; C; λ; Y]
        xp              = [up; Np; v_pretp; zp; δp]
        yp              = [θp; qp; Lp; vp; ep; Kp; Qp; pp; N_ep; ν_fp; d_fp; w_Rp; wp; L_ep; L_cp; Y_cp; Cp; λp; Yp]
        variables       = [x; y; xp; yp]
        varnames = vcat(Symbol.(x), Symbol.(y))

    # Shock
        @syms epsilon
        ex               = [epsilon]
        eta = Array([0.0; 0.0; 0.0; -σ_z; σ_δ]) # size of state space

    # Array of symbolics storing model equtions 
    f = fill(Sym("x"), 24)
    # Equilibrium conditions
        # Job creation condition -> θ
        f[1]  =  κ + K/q - β*λp/λ*(1-δbar*δ)*(w_Rp-wp-Kp+(1-s)*(κ+Kp/qp))
        # Marginal revenue product -> w_R
        f[2] = w_R - p*z*zbar/μ
        # Wage equation -> w
        f[3] = w - (ϕ*(w_R-K+θ/(1-δbar*δ)*(K+q*κ)) +(1-ϕ)*b)
        # Value of a vacancy -> Q
        f[4] = Q - (e/F)^(ξ_inv)
        # Expected discounted difference in vacancy value -> K
        f[5] = K - (Q-β*λp/λ*(1-δbar*δ)*Qp)
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
        f[12] = ν_f - β*(1-δbar*δ)*λp/λ*(ν_fp+d_fp)
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
        f[20] = v_pretp - (1-δbar*δ)*((1-q)*v+s*(1-u))
        # LOM of unemployment
        f[21] = up - ((1-(1-δbar*δ)*(θ*q))*u + (1-(1-δ*δbar)*(1-s))*(1-u))
        # LOM of firms 
        f[22] = Np - (1-δbar*δ)*(N+N_e)
        # Profits 
        #f[23] = d_f - Y_c/(N*ε)

        # Exogenous processes
        f[23]  =   log(zp) - ρ_z * log(z)
        f[24] =    log(δp) -  ρ_δ * log(δ)

    # Steady state     
        # Values 

        # Initial parameters: targets and normalizations
        N_s = 1.0
        w_s = 1.0
        z_s = 1.0
        δ_s = 1.0

        fbar = 0.41
        qbar = 0.8
        x_v = 0.20 
        labor_share = 0.66

        p_s = N_s^(1/(ε-1))
        N_es = δbar/(1-δbar)*N_s
        q_s = qbar/(1-δbar)
        θ_s = fbar/qbar
        u_s = τbar/(τbar+(1-δbar)*(θ_s*q_s))
        L_s = 1 - u_s 
        v_s = θ_s*u_s
        e_s = δbar*(v_s+1-u_s)
        v_prets = v_s - e_s 
        recruiter_share= (δbar+(ρ+δbar)*(ε-1))/(δbar+(ρ+δbar)*ε)
        w_wR = labor_share/recruiter_share
        w_Rs = w_s/(w_wR)
        surplus_ratio = (ρ+τbar)/(1-δbar)*(1/(q_s*x_v))
        K_s = (w_Rs-w_s)/(1+surplus_ratio) 
        κ   = (1-x_v)/x_v*K_s/q_s
        f_e = (μ-1)*zbar*(L_s/N_s)*(1-δbar)/(δbar*μ+ρ)
        ν_fs =p_s*f_e/μ
        d_fs = (ρ+δbar)/(1-δbar)*ν_fs
        L_es = N_es*f_e/zbar
        L_cs = L_s - L_es 
        Y_cs = p_s*zbar*L_cs
        Q_s = K_s*(1+ρ)/(ρ+δbar)
        C_s = Y_cs -  F/(1+ξ_inv)*(e_s/F)^(1+ξ_inv) -  κ*v_s*q_s 
        λ_s = C_s^(-σ)
        Y_s = Y_cs + ν_fs*N_es
        #x               = [u; N; v_pret; z] # predetermined
        #y               = [θ; q; L; v; e; K; Q; p; N_e; ν_f; d_f; w_R; w; L_e; L_c; Y_c; C; λ; Y]
        # Vector
        SS_block  = [log(x) for x in [u_s, N_s, v_prets, z_s, δ_s, θ_s, q_s, L_s, v_s, e_s, K_s, Q_s, p_s, N_es, ν_fs, d_fs, w_Rs, w_s, L_es, L_cs, Y_cs, C_s, λ_s, Y_s]]
        # vertical concatenate: represent both current and future variables
        SS = vcat(SS_block, SS_block)

        #PAR_SS = [ALPHA; BETA; DELTA; RHO; SIGMA; MUU; AA]
        PAR_SS = parameters[:]


    nx = length(x) 
    ny = length(y)
    nvar = nx + ny 
    ne = length(ex)
    nf = length(f) 

    # Procesing the model (no adjustment needed)                
        model = (parameters = parameters, estimate = estimate, estimation = position,
                npar = length(parameters), ns = length(estimate), 
                priors = priors,
                x = x, y = y, xp = xp, yp = yp, variables = variables,
                nx = nx, ny = ny, nvar = nvar,
                e = ex, eta = eta,
                ne = ne,
                f = f,
                nf = nf,
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
    ρ_δ = 0.975
    σ_δ = 0.0044

    @unpack  f_e, δ, s, z, b, ϕ, ρ, σ, ε, A, η_L, F, κ, ξ_inv = cal

    zbar = z 
    δbar = δ
    PAR     =   [f_e; s; zbar; δbar; b; ϕ; ρ; σ; ε; A; η_L; F; κ; ξ_inv; ρ_z; σ_z; ρ_δ; σ_δ ]


    # Solution functions (no adjustment needed)
        eta     =   eval_ShockVAR(PAR)
        PAR_SS  =   eval_PAR_SS(PAR)
        SS      =   eval_SS(PAR_SS)
        SS_err  =   eval_SS_error(PAR_SS, SS)
        deriv   =   eval_deriv(PAR_SS, SS)
        SS_max = maximum(abs.(SS_err))
        println("Residuals: $SS_max")

        ss = NamedTuple(zip(varnames, exp.(SS[1:nvar])))

        # @btime sol_mat = solve_model(model, deriv, eta)
        sol_mat = solve_model(model, deriv, eta)
        println("Model solved")


## Simulation
    # Impulse-Response functions
        flag_IR = true
        flag_logdev = true
        T_IR = 100
        # eta array considers all shocks simultaneously
        # For impulse responses, we can modify arrays to consider one shock at a time 
        eta_z = zero(eta)
        eta_z[4] = eta[4]

        eta_δ = zero(eta) 
        eta_δ[5] = eta[5]

        # Technology shock
        sim_IR = simulate_model(model, sol_mat, T_IR, eta_z, SS, flag_IR, flag_logdev)
        irf_df = 100 .*DataFrame(sim_IR, varnames)

        # Destruction rate shock 
        sim_IR = simulate_model(model, sol_mat, T_IR, eta_δ, SS, flag_IR, flag_logdev)
        irf_df = 100 .*DataFrame(sim_IR, varnames)

        #using Plots
        #Plots.plot(sim_IR,xlabel="Periods", ylabel= "%", yformatter=:percent)

        fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(16, 12))
        ax[1,1].plot(irf_df.u, label=:u, alpha=0.6)
        ax[1,1].plot(irf_df.v, label=:v, alpha=0.6)
        ax[1,1].plot(irf_df.θ, label=:θ, alpha=0.6)
        ax[1,1].plot(irf_df.e, label=:e, alpha=0.6)
        ax[1,1].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
        ax[1,1].legend()

        ax[1,2].plot(irf_df.C, label=:C, alpha=0.6)
        ax[1,2].plot(irf_df.Y_c, label=:Y_c, alpha=0.6)
        ax[1,2].plot(irf_df.Y, label=:Y, alpha=0.6)
        ax[1,2].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
        ax[1,2].legend()

        ax[2,1].plot(irf_df.N_e, label=:N_e, alpha=0.6)
        ax[2,1].plot(irf_df.N, label=:N, alpha=0.6)
        ax[2,1].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
        ax[2,1].legend()

        ax[2,2].plot(irf_df.z, label=:z, alpha=0.6)
        #ax[2,2].plot(irf_df.w, label=:w, alpha=0.6)
        #ax[2,2].plot(irf_df.w_R,label=:w_R, alpha=0.6)
        ax[2,2].plot(irf_df.Y, label=:Y, alpha=0.6)
        ax[2,2].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
        ax[2,2].legend()
        display(fig)



## Simulation    
        flag_IR = false
        flag_logdev = false
        T_SM = 100
        sim_SM = simulate_model(model, sol_mat, T_SM, eta, SS, flag_IR, flag_logdev)
        using Plots
        Plots.plot(sim_SM)
