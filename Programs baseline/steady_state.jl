
# =============================================================================
# Steady-State and Calibration Analysis
# =============================================================================

using Parameters, CSV, StatsBase, Statistics, Random
using NLsolve, DataFrames
using Roots, Optim, LeastSquaresOptim
using PrettyPrinting, LaTeXStrings


cd(@__DIR__)

# =============================================================================
# Matching Functions
# =============================================================================

"""
    jf(θ, A, η_L)

Job finding probability as a function of market tightness θ.
"""
jf(θ, A, η_L) = min(A * θ^(1 - η_L), 1.0)

"""
    vf(θ, A, η_L)

Vacancy filling probability as a function of market tightness θ.
"""
vf(θ, A, η_L) = min(A * θ^(-η_L), 1.0)

"""
    θ_invert(q, A, η_L)

Invert vacancy filling probability to obtain market tightness θ.
"""
θ_invert(q, A, η_L) = (q / A)^(-1 / η_L)

# =============================================================================
# Parameter Struct
# =============================================================================

@with_kw mutable struct ParaCalib
    f_e::Float64 = 1.0                 # Entry cost
    τ::Float64 = 0.0305                # Worker separation rate
    δ::Float64 = 0.005                 # Product destruction rate
    z::Float64 = 1.0                   # Technology level
    b::Float64 = 0.71                  # Unemployment insurance
    ϕ::Float64 = 0.5                   # Bargaining power
    r::Float64 = 0.04 / 12             # Rate of time preference (monthly)
    σ::Float64 = 1.0                   # Inverse intertemporal elasticity of substitution
    ε::Float64 = 4.0                   # Elasticity of substitution
    A::Float64 = 0.5631                # Job matching level parameter
    η_L::Float64 = 0.6                 # Elasticity of matching function w.r.t. unemployment
    F::Float64 = 0.0007                # Mass of recruiters
    κ::Float64 = 0.2                   # Fixed matching cost
    ξ_inv::Float64 = 1.0               # Inverse elasticity of entry to vacancy value
    x_m::Float64 = 113.16              # Upper bound on cost distribution

    # New parameters related to firm heterogeneity
    ψ::Float64 = 1.5                   # Curvature related to power law of continuation cost
    f_m::Float64 = 10.0                 # Max value of cost draw              

    # Derived parameter: worker separation rate for steady-state
    #s = (τ - δ) / (1 - δ)
end

# =============================================================================
# Equilibrium Functions
# =============================================================================

"""
    F(x)
Continuation cost-draw cdf given maximal cost
"""
function F(x, f_m, ψ)
    (a_m <=0 || ψ<=0) && throw(ArgumentError("a_m, k must be > 0 "))
    if x < f_m
        out = (x/f_m)^ψ
    else
        out = 1.0
    end 
    return out
end 

"""
F_inv(x_m, f_m, ψ)
Inverse cdf: Given probability F, finds associated value x
"""
function F_inv(F, f_m, ψ)
    # Solve F = (x/f_m)^ψ for x
    (f_m <=0 || ψ<=0) && throw(ArgumentError("f_m, ψ must be > 0 "))
    x = F^(1/ψ)*f_m
    return x 
end 

"""
    L_fun(θ, para)
Steady-state employment as a function of market tightness θ, δ_e and parameters.
"""
function L_fun(θ, δ_e, para)
    @unpack δ, s, A, η_L = para
    f = jf(θ, A, η_L)
    τ = 1 - (1 - δ_e) * (1 - s)
    return (1 - δ_e) * f / (τ + (1 - δ_e) * f)
end

"""
    e_fun(θ, para)
New vacancy rate as a function of market tightness θ and parameters.
"""
function e_fun(θ, δ_e, para)
    @unpack δ, A, η_L = para
    f = jf(θ, A, η_L)
    τ = 1 - (1 - δ_e) * (1 - s)
    return δ_e * (θ * τ + (1 - δ_e) * f) / (τ + (1 - δ_e) * f)
end

"""
    K_fun(θ, para)
Vacancy value as a function of market tightness θ and parameters.
"""
function K_fun(θ, δ_e, para)
    @unpack ρ, δ, x_m, ξ_inv = para
    e = e_fun(θ, δ_e, para)
    # Value of vacancy
    Q = (e)^ξ_inv * x_m
    return (r + δ_e) / (1 + r) * Q
end

"""
    w_fun(θ, N, para)

Wage function as a function of market tightness θ, number of businesses N, and parameters.
"""
function w_fun(θ, N, δ_e, a_tilde, para)
    @unpack ϕ, b, z, ε, δ, A, η_L, κ = para
    q = vf(θ, A, η_L)
    μ = ε / (ε - 1)
    w_int = N^(1 / (ε - 1)) * z*a_tilde/ μ
    K = K_fun(θ, δ_e,  para)
    return (1 - ϕ) * b + ϕ * (w_int - K + θ * (K + q * κ))
end


# =============================================================================
# Steady-State Solver
# =============================================================================

"""
    steady_state(para; init=0.51)

Compute steady-state statistics given parameters.
Returns a NamedTuple of all key statistics.
"""
function steady_state(para; init=0.51)
    @unpack f_e, τ, δ, s, z, b, ϕ, ρ, σ, ε, A, η_L, κ, ξ_inv, f_m, ψ, s = para
    μ = ε / (ε - 1)
    #θ = θ_fun(para, init_value=init)
    out = zeros(2)

    function loss(x)
        θ, end_dest = abs.(x) #ensure non-negativity
        # Aggregate destruction rate
        δ_e = δ + end_dest - δ*end_dest   #(1-δ_e)=(1-δ)*(1-end_dest)
        f = jf(θ, A, η_L)
        q = vf(θ, A, η_L)
        K = K_fun(θ, δ_e, para)
        L = L_fun(θ, δ_e, para)
        u = 1 - L
        lhs_jcc = (κ + K / q) * (r + τ + (1 - δ_e) * ϕ * q * θ)
        N = (μ - 1) * z * L * (1 - δ_e) / (f_e * (δ_e * μ + r))

        ρ = N^(1 / (ε - 1))
        w_int = ρ * z / μ
        rhs_jcc = (1 - δ_e) * (1 - ϕ) * (w_int - K - b)

        N_e = δ_e*N/(1-δ_e)
        L_c = (r + δ_e) * L / (δ_e * μ + r)
        # Consumption output
        Y_c = ρ * z * L_c
        # Cutoff
        x_c = Y_c/(ε*N) + ρ*f_e/μ
        surv_prob = F(x_c, f_m, ψ) #∈ (0, 1)
        δ_e_new = 1.0 - surv_prob*(1-δ) #(1-δ_e_new) = surv_prob*(1-δ)

        vars = (; δ_e, f, q, θ, u, ρ, w_int, K, N, N_e, L_c, Y_c, x_c)

        out[1] = (rhs_jcc-lhs_jcc)/(rhs_jcc+lhs_jcc)
        out[2] = δ_e - δ_e_new
        return out, vars 
    end 

    sol = LeastSquaresOptim.optimize(x -> loss(x)[1], [0.51, δ/2.0], Dogleg())
    println("converged=$(sol.converged) at root=$(sol.minimizer) in " *
        "$(sol.iterations) iterations and $(sol.f_calls) function calls")
    
    var = loss(sol.minimizer)[2]
    @unpack δ_e, f, q, θ, u, ρ, w_int, K, N, N_e, L_c, Y_c, x_c = var
    L_e = L-L_c

    # Labor market variables
    v = θ * u
    e = δ_e * (v + 1 - u)
    Q = e^ξ_inv * x_m
    X_v = e * (1 / (1 + ξ_inv)) * Q
    v_pret = v - e

    # Relative price, businesses, and values
    ν_f = ρ * f_e / μ
    d_f = (ρ + δ_e) / (1 - δ_e) * ν_f

    # Wages
    w = ϕ * (w_int - K + θ * (K + q * κ)) + (1 - ϕ) * b
    # Sectoral labor 
    L_e = δ_e * (μ - 1) * L / (δ * μ + r)
    # Total recruiting costs
    X = X_v + κ * v * q

    # Total (stochastic) fixed costs 
    X_c = N*(ψ/(ψ+1))*x_c 

    # Output and shares
    C = Y_c - X - X_c
    Y = C + ν_f * N_e
    J = Q + (1 + r) / (1 - δ) * K / q
    M = Q * v + J * L + (N + N_e) * ν_f
    labor_prod = Y / (ρ*L)

    labor_share = w * L / Y
    cons_share = C / Y
    inv_new_firm_share = ν_f * N_e / Y
    vacancy_share = X / Y
    fixed_cost_share = X_c/Y
    sunk_vac_cost_share = X_v / Y
    entrant_share = e / v

    x_v = (K / q) / (κ + K / q)
    search_wedge = w / w_int
    recruiter_share = w_int * L / Y
    ann_int_rate = (1 + r)^12 - 1

    return (
        θ=θ, N=N, f=f, q=q, u=u, v=v, v_pret=v_pret, e=e, K=K, ρ=ρ, N_e=N_e,
        ν_f=ν_f, d_f=d_f, w_int=w_int, w=w, L=L, L_e=L_e, L_c=L_c, Q=Q, J=J,
        X_v=X_v, X=X, C=C, Y_c=Y_c, Y=Y, labor_share=labor_share,
        labor_prod=labor_prod, cons_share=cons_share, inv_new_firm_share=inv_new_firm_share,
        vacancy_share=vacancy_share, sunk_vac_cost_share=sunk_vac_cost_share, M=M,
        entrant_share=entrant_share, x_v=x_v, search_wedge=search_wedge,
        recruiter_share=recruiter_share, μ=μ, ann_int_rate=ann_int_rate
    )
end

# =============================================================================
# Calibration
# =============================================================================

# Default calibration targets
targets = (
    X_Y=0.015, dest_ann=0.0754, dest_end_frac=0.5, f=0.41, η_L=0.6, q=0.8, sep=0.031, b_ratio=0.71,
    x_v=1.0, ξ_inv=1, r_ann=0.04, ε=4.3, σ=1.0, N=1.0, w=1.0, f_m=1.0, ψ=1.5)

"""
    calibrate_shares(targets)

Calibrate model parameters to match empirical targets.
Returns a NamedTuple of calibrated parameters.
"""
function calibrate_shares(targets)
    @unpack X_Y, dest_ann, dest_end_frac, f, η_L, q, sep, b_ratio, x_v, ξ_inv, ε, r_ann, σ, N, w, f_m, ψ = targets

    # Aggregate monthly destruction rate
    δ_e = 1 - (1 - dest_ann)^(1 / 12)
    # Exogenous destruction rate
    δ = (1.0-dest_end_frac)*δ_e 
    # Probability of surviving fixed cost shock
    surv_prob = (1-δ_e)/(1-δ)

    μ = ε/(ε-1)

    # Purely endogenous destruction threshold

    τ = sep
    s = (τ - δ_e) / (1 - δ_e)
    r = (1 + r_ann)^(1 / 12) - 1

    # Correct job finding and vacancy filling probabilities
    f = f / (1 - δ_e)
    q = q / (1 - δ_e)

    θ = f / q
    u = τ / (τ + (1 - δ_e) * f)
    v = θ * u
    L = 1 - u

    # Matching function level parameter
    A = f / θ^(1 - η_L)

    e = δ_e * (v + 1 - u)
    N_e = δ_e / (1 - δ_e) * N
    b = b_ratio * w

    # Sectoral labor 
    L_c = (r+δ_e)*L/(δ_e*μ+r)
    L_e = L - L_c 
    #@assert L_e == δ_e*(μ-1)*L/(δ_e*μ+r)
    #Retailer profit share (would be 1/ε in absence of fixed costs)
    #π_s = (1/ε - f_r) #Xc_Yc



    surplus_ratio = (r + τ) / (1 - δ_e) * (1 / (q * x_v))
    ρ = N^(1 / (ε - 1))
    μ = ε / (ε - 1)
    #recruiter_share = (δ + (ρ + δ) * (ε - 1)) / (δ + (ρ + δ) * ε)
    #recruiter_share = ((1-π_s)*ρ+δ_e)/(ρ+δ_e*(1+π_s))

    X_c_share = 0.20

    function loss(Q)
        Q = abs(Q)
        K = Q * (r + δ_e) / (1 + r)
        κ = (1 - x_v) / x_v * K / q
        # Aggregate recruiting costs: (1) sunk and (2) fixed matching costs
        X = e / (1 + ξ_inv) * Q + κ * q * v
        w_int = surplus_ratio * K + w + K
        ϕ = (w - b) / (w_int - K + θ * (K + q * κ) - b)

        z = (μ / ρ) * w_int # ρ = μ*w_int/z
        f_e = (μ - 1) * z * (L / N) * (1 - δ_e) / (δ_e * μ + r) # rearrange resource constraint curve
        ν_f = ρ * f_e / μ
        Y_c = ρ * z * L_c
        # Cutoff
        x_c = Y_c/(ε*N) +ν_f 
        X_c = X_c_share*Y_c
        # Aggregate fixed costs paid 
        #X_c = N*ψ/(ψ+1)*x_c
        # Consumption and output
        C = Y_c - X - X_c # Net out intermediate goods X and X_c
        Y_new = C + ν_f * N_e
        Y = 1 / X_Y * X
        return (Y - Y_new) / (Y + Y_new), (;w_int, κ, z, f_e, K, x_c, X_c, X, C, Y_c, Y)
    end



    Q = fzero(x -> loss(x)[1], 0.1)
    Q = abs(Q)
    out = loss(Q)[2]
    @unpack w_int, κ, z, f_e, K, x_c, X_c, X, C, Y_c, Y = out
    X_Y = abs(X_Y)

    # Solve for ψ 
    ψ_coeff = X_c/(N*x_c)
    ψ = ψ_coeff/(1-ψ_coeff)
    # Destruction rate 
    # Implied power law parameter: (x_c/f_m)^ψ = surv_prob
    f_m = x_c/surv_prob^(1/ψ)
    @assert (F(x_c, f_m, ψ) - surv_prob) == 0.0


    ϕ = (w - b) / (w_int - K + θ * (K + q * κ) - b)
    x_m = Q / e^ξ_inv

    return (
        f_e=f_e, τ=τ, δ=δ, z=z, b=b, ϕ=ϕ, ρ=ρ, σ=σ, ε=ε, A=A, η_L=η_L,
        κ=κ, ξ_inv=ξ_inv, x_m=x_m, s=s, ψ=ψ, f_m=f_m
    )
end

# =============================================================================
# Curves: Job Creation & Resource Constraint
# =============================================================================

"""
    N_jcc(θ, para)

Job creation curve: N as a function of θ.
"""
function N_jcc(θ, para)
    @unpack ρ, τ, A, η_L, δ, ϕ, z, ε, b, κ = para
    q = vf(θ, A, η_L)
    K = K_fun(θ, para)
    μ = ε / (ε - 1)
    w_int = (1 / ((1 - δ) * (1 - ϕ))) * (κ * q + K) * ((ρ + τ) / q + (1 - δ) * ϕ * θ) + K + b
    p = w_int * (μ / z)
    N = p^(ε - 1)
    return N
end

"""
    N_res(θ, para)

Resource constraint curve: N as a function of θ.
"""
function N_res(θ, para)
    @unpack ρ, δ, τ, A, η_L, f_e, ε, z = para
    μ = ε / (ε - 1)
    f = jf(θ, A, η_L)
    u = τ / (τ + (1 - δ) * f)
    L = 1 - u
    N = (μ - 1) * z * L * (1 - δ) / (f_e * (δ * μ + ρ))
    return N
end

