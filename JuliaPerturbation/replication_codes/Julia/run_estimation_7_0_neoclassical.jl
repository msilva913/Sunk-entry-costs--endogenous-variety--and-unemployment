# (c) Alvaro Salazar-Perez and Hernán D. Seoane
# "Perturbating and Estimating DSGE models in Julia
# This version 2023

# If you use these codes, please cite us.
# These codes are free and without any guarantee. Use them at your own risk.
# If you find typos, please let us know.
cd(@__DIR__)
#v1.7+
using MKL

#v1.7- 
#BLAS.vendor() 
#:mkl

include("solution_functions_7_0.jl")

## Model
        # Adjustment
                flag_order      = 1 # perturbation ordr
                flag_deviation  = true 
                flag_SSsolver   = false

        # Parameters
                @vars α β δ ρ σ μ A
                parameters      = [α; β; δ; ρ; σ; μ; A]
                estimate        = [ρ; μ]
                position        = [4;6]
                # (initial value, lower bound, upper bound, distribution)
                priors          = (prior_rho = (iv = .9 , lb = 0.0 , ub = 1.0, d = truncated(Normal(0.75, 0.25))),
                                   prior_muu = (iv = .05, lb = 0.0, ub = 1e6, d = InverseGamma(5.0, 0.25)))

        # Variables
                @vars k kp a ap c cp n np yy yyp r rp ii iip
                x               = [k; a] # states at t
                y               = [c; n; r; yy; ii] # jumpers at t
                xp              = [kp; ap] # states at t+1
                yp              = [cp; np; rp; yyp; iip] #jumpers at t+1
                variables       = [x; y; xp; yp] #group variables

        # Shock
                @vars epsilon
                e               = [epsilon] # vector of exogenous shocks
                eta             = Array([0.0; μ])
                
        # Equilibrium conditions
                f1  =   c + kp - (1-δ) * k - a * k^α * n^(1-α)
                f2  =   c^(-σ) - β * cp^(-σ) * (ap *α* kp^(α-1) * np^(1-α) + 1 - δ)
                f3  =   A - c^(-σ) * a * (1-α) * k^α * n^(-α) 
                f4  =   log(ap) - ρ * log(a)
                f5  =   r - a * α * k^(α-1) * n^(1-α)
                f6  =   yy - a * k^α * n^(1-α)
                f7  =   ii - (kp - (1-δ) * k)

                f   =   [f1;f2;f3;f4;f5;f6;f7]

        # Steady State
                # Values 
                A   =   1.0
                N   =   2/3
                K   =   (α/(1/β-1+δ))^(1/(1-α))*N
                C   =   A * K^(α) * N^(1-α) - δ*K
                R   =   A * α * K^(α-1.0) * N^(1-α)
                YY  =   A * K^α * N^(1-α)
                II  =   δ*K
                #Vector
                SS = [
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
                AA = C^(-σ) * (1-α)*A*K^(α)*N^(-α)
                PAR_SS = [α; β; δ; ρ; σ; μ; A]
        
        # Procesing the model (no adjustment needed) 
                model = (parameters = parameters, #total parameters
                 estimate = estimate, #parameters to be estimated
                 estimation = position, #position of estimated parameters
                 npar = length(parameters), ns = length(estimate), #length of parameters and estimates
                  priors = priors, #prior distributions
                  x = x, y = y, xp = xp, yp = yp, variables = variables,
                        nx = length(x), # length of states
                         ny = length(y), # length of jumpers
                          nvar = length(x) + length(y), #total length of vars
                        e = e, eta = eta,
                        ne = length(e),
                        f = f,
                        nf = length(f),
                        SS = SS, PAR_SS = PAR_SS,
                        flag_order = flag_order, flag_deviation = flag_deviation, flag_SSsolver = flag_SSsolver)

                process_model(model)

## Solution
        # Parametrization
                α  =   0.30
                β   =   0.95
                δ  =   1.00
                ρ   =   0.90
                σ  =   2.00
                μ   =   0.05

                PAR     =   [α; β; δ; ρ; σ; μ]
        
        # Solution functions (no adjustment needed)
                eta     =   eval_ShockVAR(PAR)
                PAR_SS  =   eval_PAR_SS(PAR)
                SS      =   eval_SS(PAR_SS)
                deriv   =   eval_deriv(PAR_SS, SS)

                # y_t = g(x_t, σ), solution for jump varibles 
                #x_{t+1} = h(x_t, σ) + ησϵ_{t+1}, solution for state variables

                @btime sol_mat = solve_model(model, deriv, eta) # Solution g and h matrices
                sol_mat = solve_model(model, deriv, eta)

## Estimation
        # Simulate a sample for estimation                
                nsim = 200
                # nsim = 500
                discard = 10000
                flag_IR = false
                flag_logdev = true
                simulation_logdev = simulate_model(model, sol_mat, discard+nsim, eta, SS, flag_IR, flag_logdev)
                data = simulation_logdev[discard+1:discard+nsim,model.nx+1:model.nx+3]
        
        # Estimation
                c       = 0.1           # step size
                npart   = 500           # of particles
                # npart   = 2^13        # of particles
                nphi    = 100           # of stage
                lam     = 3             # bending coeff
                #@time estimation_results = smc_rwmh(model, data, PAR, PAR_SS, eta, SS, deriv, sol_mat, c, npart, nphi, lam)