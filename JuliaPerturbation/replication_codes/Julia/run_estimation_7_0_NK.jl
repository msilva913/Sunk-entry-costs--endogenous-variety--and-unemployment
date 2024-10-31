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
    # Adjustment
        flag_order      = 1
        flag_deviation  = true
        flag_SSsolver   = false

    # Parameters
    "
    α: Capital share in production 
    β: Discount factor
    γ: Risk aversion parameter (CRRA)
    δ: Capital depreciation 
    ρ_z: TFP persistence coefficient
    σ_z: TFP standard deviation 
    Z: TFP mean 
    η: Coefficient of elasticity of subs. b/w varieties
    θ: Calvo probability of being able to switch prices 
    ϕ_Π: Interest rate response to inflation 
    ϕ_Y: Interest rate response to output 
    ϕ_R: Smoothing coefficient of MP rule 
    "
        # List all parameters in model to generate them as symbol
        @vars θ β γ δ η α ξ ρ ϕ_R ϕ_Π ϕ_Y Π_star YSS RRS σ_a
        parameters  =   [θ; β; γ; δ; η; α; ξ; ρ; ϕ_R; ϕ_Π; ϕ_Y; Π_star; YSS; RRS; σ_a]
        # List all parameters that need to be estimatd
        estimate    =   [ρ; σ_a]
        #position = [8; 15]
        position = [findfirst(==(param), parameters) for param in estimate]
        # Specify priors: initial value, lower bound, upper bound, distribution
        priors = (prior_rho = (iv = .9 , lb = 0.0 , ub = 1.0, d = truncated(Normal(0.75, 0.25),0.0,1.0)),
                  prior_sigma = (iv = .05, lb = 0.0, ub = 1e6, d = InverseGamma(5.0, 0.25)))

    # Variables
    """
        c, c: consumption
        h, hp: habit formation 
        """
        # List all variables in model to generate them as symbols
        @vars c cp h hp l lp w wp mc mcp s sp ph php ppi ppip rt rtp x1 x1p x2 x2p br brp brb brbp z zp k kp ii iip u up
        x = [s; z; brb; k] #states
        y = [c; h; ppi; rt; br; w; l; mc; ph; x1; x2; ii; u] #jump variables
        xp = [sp; zp; brbp; kp] # next-period states
        yp = [cp; hp; ppip; rtp; brp; wp; lp; mcp; php; x1p; x2p; iip; up] #next-period jumpers
        variables = [x; y; xp; yp]

    # Shock
        # List all exogenous shocks to generate them as symbols 
        @vars epsilon
        e = [epsilon]
        # Standard deviations
        eta = Array([0.0; σ_a; 0.0; 0.0])

    # Equilibrium conditions as a symbolic function 
        f1 = -kp + (1-δ)*k + ii
        f2  = l - (c*(1-h)^γ)^(-ξ) * (1-h)^γ
        f3  = β*lp/ppip - l*rt
        f4  = 1/br - rt
        f5  = -sp + (1-θ)*(ph^(-η)) + θ*(ppi^η)*s
        f6  = z*h^(1-α)*k^(α) - sp*(c+ii)
        f7  = -1 + (1-θ)*(ph^(1-η)) + θ*(ppi^(η-1))
        f8  = -x1 + (ph^(1-η))*c*((η-1)/η) + θ*rt*((ph/php)^(1-η))*(ppip^(η))*x1p
        f9  = -x2 + (ph^(-η))*c*mc + θ*rt*((ph/php)^(-η))*(ppip^(η+1))*x2p
        f10 = x2-x1
        f11 = w - mc*z*(1-α)*h^(-α)*k^(α)
        f12 = w - ((c*(1-h)^γ)^(-ξ)*γ*(1-h)^(γ-1)*c)/((c*(1-h)^γ)^(-ξ)*(1-h)^γ)
        f13a = u - mc*z*α*h^(1-α)*k^(α-1)
        f13b = l - β*(up+1-δ)*lp
        f14a = log(br/RRS) - ϕ_R*log(brb/RRS) - ϕ_Π*log(ppi/Π_star) - ϕ_Y*log((sp*(c+ii))/YSS)
        f14b = brbp - br
        f15 = log(zp) - ρ * log(z)

        f = [f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12;f13a;f13b;f14a;f14b;f15]

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
        Z           =       1
        R           =       1/RRS
        RR          =       RRS
        PH          =       ((1-θ*Π_star^(η-1))/(1-θ))^(1/(1-η))
        PPI         =       Π_star
        S           =       (1-θ)*PH^(-η)/(1-θ*Π_star^(η))
        MC          =       ((1-θ*RRS^(-1)*Π_star^(η+1))*((η-1)/η)*PH)/(1-θ*RRS^(-1)*Π_star^(η))
        U           =       1/β - 1 + δ
        H           =       ((1-α)*MC*Z*((α*Z*MC)/U)^(α/(1-α)))/((1-α)*MC*Z*((α*Z*MC)/U)^(α/(1-α)) + (γ*Z*((α*Z*MC)/U)^(α/(1-α)) - S*δ*((α*Z*MC)/U)^(1/(1-α)))/S);
        K           =       ((α*Z*MC)/U)^(1/(1-α))*H
        W           =       (1-α)*MC*Z*(K/H)^α
        Y           =       Z*H^(1-α)*K^α
        YSS         =       Y
        C           =       (Y - S*δ*K)/S
        II          =       δ*K
        X1          =       (PH^(1-η)*Y*S^(-1)*((η-1)/η))/(1-θ*RRS^(-1)*Π_star^η)
        X2          =       X1
        LAM         =       C^(-ξ)*(1-H)^(γ*(1-ξ))
        # Vector
        SS = [#arrange in order
                log(S);
                log(Z);
                log(RR);
                log(K);
                log(C);
                log(H);
                log(PPI);
                log(R);
                log(RR);
                log(W);
                log(LAM);
                log(MC);
                log(PH);
                log(X1);
                log(X2);
                log(II);
                log(U);
                log(S);
                log(Z);
                log(RR);
                log(K);
                log(C);
                log(H);
                log(PPI);
                log(R);
                log(RR);
                log(W);
                log(LAM);
                log(MC);
                log(PH);
                log(X1);
                log(X2);
                log(II);
                log(U);
            ]

        #PAR_SS = [THETA; BETTA; GAMMA; DELTA; ETA; ALPHA; XI; RHO; PHI_R; PHI_PPI; PHI_Y; PPS; YSS; RRS; SIGMAA]
        PAR_SS = parameters[:]

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
    # Parametrization #
    #@vars θ β γ δ η α ξ ρ ϕ_R ϕ_Π ϕ_Y PPS YSS RRS σ_a
        θ      =       0.80
        β      =       1.04^(-1/4)
        γ      =       3.6133
        δ      =       0.01
        η        =       5.00
        α      =       0.30
        ξ         =       2.00
        ρ        =       0.8556
        Z          =       1.00
        ϕ_Π    =       3.00
        ϕ_Y      =       0.01
        ϕ_R      =       0.80
        Π_star        =       1.042^(1/4)
        RRS        =       Π_star/β
        σ_a     =       0.05

        Z          =       1.00
        R          =       1.00/RRS
        RR         =       RRS
        PH         =       ((1.00-θ*Π_star^(η-1.00))/(1.00-θ))^(1.00/(1.00-η))
        PPI        =       Π_star
        S          =       (1.00-θ)*PH^(-η)/(1.00-θ*Π_star^(η))
        MC         =       ((1.00-θ*RRS^(-1.00)*Π_star^(η+1.00))*((η-1.00)/η)*PH)/(1.00-θ*RRS^(-1.00)*Π_star^(η))
        U          =       1.00/β - 1.00 + δ
        H           =       ((1-α)*MC*Z*((α*Z*MC)/U)^(α/(1-α)))/((1-α)*MC*Z*((α*Z*MC)/U)^(α/(1-α)) + (γ*Z*((α*Z*MC)/U)^(α/(1-α)) - S*δ*((α*Z*MC)/U)^(1/(1-α)))/S)
        K           =       ((α*Z*MC)/U)^(1/(1-α))*H
        W          =       (1.00-α)*MC*Z*(K/H)^α
        YSS        =       Z*H^(1.00-α)*K^α

        PAR = [θ; β; γ; δ; η; α; ξ; ρ; ϕ_R; ϕ_Π; ϕ_Y; Π_star; YSS; RRS; σ_a]

    # Solution functions (no adjustment needed)    
        eta     =   eval_ShockVAR(PAR)
        PAR_SS  =   eval_PAR_SS(PAR)
        SS      =   eval_SS(PAR_SS)
        deriv   =   eval_deriv(PAR_SS, SS)
        
        # @btime sol_mat = solve_model(model, deriv, eta)
        sol_mat = solve_model(model, deriv, eta)

## Estimation
    # Simulate a sample for estimation         
        nsim = 200
        # nsim = 500
        discard = 10000
        flag_IR = false
        flag_logdev = true
        simulation_logdev = simulate_model(model, sol_mat, discard+nsim, eta, SS, flag_IR, flag_logdev)
        # Extract c, h, infl
        data = simulation_logdev[discard+1:discard+nsim,model.nx+1:model.nx+3]
    
        irf = simulate_model(model, sol_mat, 100, eta, SS, true, flag_logdev)
    # Estimation
        c       = 0.1
        #npart   = 500           # of particles
        npart = 1000
        # npart   = 2^13        # of particles
        nphi    = 100           # of stage
        lam     = 3             # bending coeff
        #@time estimation_results = smc_rwmh(model, data, PAR, PAR_SS, eta, SS, deriv, sol_mat, c, npart, nphi, lam)
        para, wght, mu, sig = estimation_results = smc_rwmh(model, data, PAR, PAR_SS, eta, SS, deriv, sol_mat, c, npart, nphi, lam)

    # Mode
    indices = argmax(wght[:,1:2], dims=1)
    mode = para(indices)

    using PyPlot
    fig, ax = plt.subplots(ncols=2)
    perm = sortperm(para[:,1])
    ax[1].plot(para[:,1][perm], wght[:,1][perm], label="Density of ρ")
    perm = sortperm(para[:,2])
    ax[2].plot(para[:,2][perm], wght[:,2][perm], label="Density of σ")
    display(fig)
    plt.show()