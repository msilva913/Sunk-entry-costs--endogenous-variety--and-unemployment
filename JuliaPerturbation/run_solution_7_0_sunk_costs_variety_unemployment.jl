# (c) Alvaro Salazar-Perez and Hernán D. Seoane
# "Perturbating and Estimating DSGE models in Julia
# This version 2023

# If you use these codes, please cite us.
# These codes are free and without any guarantee. Use them at your own risk.
# If you find typos, please let us know.

#v1.7+
using MKL

#v1.7- 
#BLAS.vendor() 
#:mkl

include("solution_functions_7_0.jl")


## Model
    # Adjustments
        flag_order      = 1
        flag_deviation  = true
        flag_SSsolver   = false

    # Parameters
        @vars f_e τ δ s zbar b ϕ ρ σ ε A η_L F κ ξ_inv ρ_z μ_z
        parameters      = [f_e; τ; δ; s; z; b; ϕ; ρ; σ; ε; A; η_L; F; κ; ξ_inv; ρ_z; μ_z]
        estimate        = []
        position        = []
        priors          = (;)

        # Transformations
        β = 1/(1+ρ)
        ξ = 1/ξ_inv
        τ = 1 - (1-δ)*(1-s)
        μ = ε/(ε-1)

    # Variables
        @vars  z u θ f q K L u v v_pret e K Q X_v N p N_e ν_f d_f w_R w L_e L_c Y_c X C λ Y labor_share
        @vars zp θp fp qp Kp Lp up vp v_pretp ep Kp Qp X_vp Np pp N_ep ν_fp d_fp w_Rp wp L_ep L_cp Y_cp Xp Cp λp Yp labor_sharep
       
        
        x               = [u; N; v_pret; z] # predetermined
        y               = [θ; f; q; K; L; v; e; K; Q; X_v; N; p; N_e; ν_f; d_f; w_R; w; L_e; L_c; X; C; λ; Y; labor_share]
        xp              = [up; Np; v_pretp; zp]
        yp              = [θp; fp; qp; Kp; Lp; vp; ep; Kp; Qp; X_vp; Np; pp; N_ep; ν_fp; d_fp; w_Rp; wp; L_ep; L_cp; Xp; Cp; λp; Yp; labor_sharep]
        variables       = [x; y; xp; yp]

    # Shock
        @vars epsilon
        e               = [epsilon]
        eta             = Array([0.0; μ_z])

    # Equilibrium conditions
        # Job creation condition
        f1  =  κ + K/q - β*λp/λ*(1-δ)*(w_Rp-wp-Kp+(1-s)*(κ+Kp/qp))
        # Marginal revenue product
        f2 = w_R - p*z*zbar/μ
        # Wage equation 
        f3 = w - (ϕ*(w_R-K+θ/(1-δ)*(K+q*κ)) +(1-ϕ)*b)
        # Value of a vacancy
        f4 = Q - (e/F)^(ξ_inv)
        # Expected discounted difference in vacancy value 
        f5 = K - (Q-β*λp/λ*(1-δ)*Qp)
        # Market tightness 
        f6 = θ - v/u 
        # Job finding probability 
        f7 = f -A*θ^(1-η_L) 
        # Vacancy filling probability 
        f8 = q - A*θ^(-η_L) 
        # Aggregate labor 
        f9 = L - (1-u)
        # Composition of labor 
        f10 = L - (L_c+L_e) 
        # Relative price 
        f11 = p - N^(1/(ε-1))
        # Lagrangian multiplier 
        f12 = λ - C^(-σ)
        # Firm value 
        f13 = ν_f - β*(1-δ)*λp/λ*(ν_fp+d_fp)
        # Retail output: resources 
        f14 = Y_c - p*z*zbar*L_c
        # Retail output: expenditure 
        f15 = Y_c - (C+F/(1+ξ_inv)*(e/F)^(1+ξ_inv)+κ*v*q)
        # Business entrants
        f17 = N_e - L_e*z*zbar/f_e 
        # Firm value relative to price 
        f18 = ν_f - p*f_e/μ 
        # Output = expenditure
        f19 = Y - (Y_c+ν_f*N_e)
        # Output = income 
        f20 = Y - (w_R*L+N*d_f)
        # LOM of vacancies 
        f21 = v - (v_pret + e)
        # Predetermined vacancies 
        f22 = v_pretp - (1-δ)*((1-q)*v+s*(1-u))
        # LOM of unemployment
        f23 = up - ((1-(1-δ)*f)*u + τ*(1-u))
        # LOM of firms 
        f24 = Np - (1-δ)*(N+N_e)
        # Labor share 
        f25 = labor_share - w*L/Y 

        # Exogenous processes
        f26  =   log(zp) - ρ_z * log(z)

        f   =   [f1;f2;f3;f4;f5;f6;f7]

    # Steady state     
        # Values 
            A   =   1.0
            N   =   2/3
            K   =   (ALPHA/(1/BETA-1+DELTA))^(1/(1-ALPHA))*N
            C   =   A * K^(ALPHA) * N^(1-ALPHA) - DELTA*K
            R   =   A * ALPHA * K^(ALPHA-1.0) * N^(1-ALPHA)
            YY  =   A * K^ALPHA * N^(1-ALPHA)
            II  =   DELTA*K
        # Vector
            SS =    [
                    log(K); # k
                    log(A); # a
                    log(C); # c
                    log(N);
                    log(R);
                    log(YY);
                    log(II);
                    log(K); # kp
                    log(A); # ap
                    log(C); # cp
                    log(N);
                    log(R);
                    log(YY);
                    log(II)
                    ]

        # Parameters to adjust
            AA = C^(-SIGMA) * (1-ALPHA)*A*K^(ALPHA)*N^(-ALPHA)
            PAR_SS = [ALPHA; BETA; DELTA; RHO; SIGMA; MUU; AA]

    # Procesing the model (no adjustment needed)                
        model = (parameters = parameters, estimate = estimate, estimation = position,
                npar = length(parameters), ns = length(estimate), priors = priors,
                x = x, y = y, xp = xp, yp = yp, variables = variables,
                nx = length(x), ny = length(y), nvar = length(x) + length(y),
                e = e, eta = eta,
                ne = length(e),
                f = f,
                nf = length(f),
                SS = SS, PAR_SS = PAR_SS,
                flag_order = flag_order, flag_deviation = flag_deviation, flag_SSsolver = flag_SSsolver)
        process_model(model)


## Solution
    # Parametrization
        ALPHA  =   0.30
        BETA   =   0.95
        DELTA  =   1.00
        RHO    =   0.90
        SIGMA  =   2.00
        MUU    =   0.05

        PAR     =   [ALPHA; BETA; DELTA; RHO; SIGMA; MUU]


    # Solution functions (no adjustment needed)
        eta     =   eval_ShockVAR(PAR)
        PAR_SS  =   eval_PAR_SS(PAR)
        SS      =   eval_SS(PAR_SS)
        SS_err  =   eval_SS_error(PAR_SS, SS)
        deriv   =   eval_deriv(PAR_SS, SS)
        println("Residuals: $SS_err")

        # @btime sol_mat = solve_model(model, deriv, eta)
        sol_mat = solve_model(model, deriv, eta)
        println("Model solved")


## Simulation
    # Impulse-Response functions
        flag_IR = true
        flag_logdev = true
        T_IR = 30
        sim_IR = simulate_model(model, sol_mat, T_IR, eta, SS, flag_IR, flag_logdev)
        using Plots
        plot(sim_IR,xlabel="Periods", ylabel= "%", yformatter=:percent)


## Simulation    
        flag_IR = false
        flag_logdev = false
        T_SM = 100
        sim_SM = simulate_model(model, sol_mat, T_SM, eta, SS, flag_IR, flag_logdev)
        using Plots
        plot(sim_SM)
