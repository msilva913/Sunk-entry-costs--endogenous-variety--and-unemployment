
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
    ρ::Float64 = 0.04 / 12             # Rate of time preference (monthly)
    σ::Float64 = 1.0                   # Inverse intertemporal elasticity of substitution
    ε::Float64 = 4.0                   # Elasticity of substitution
    A::Float64 = 0.5631                # Job matching level parameter
    η_L::Float64 = 0.6                 # Elasticity of matching function w.r.t. unemployment
    F::Float64 = 0.0007                # Mass of recruiters
    κ::Float64 = 0.2                   # Fixed matching cost
    ξ_inv::Float64 = 1.0               # Inverse elasticity of entry to vacancy value
    x_m::Float64 = 113.16              # Upper bound on cost distribution

    # New parameters related to firm heterogeneity
    f::Float64 = 2.0                   # Fixed operating cost 
    a_m::Float64 = 1.0                 # Tail parameter of Pareto distribution => minimum productivity draw
    k::Float64 = 3.4                 

    #Notes:
    #1) We must have k >  ε-1 to have a well-define price index 
    #2) a_m=1.0 just fixes units, from which we can obtain a_star and a_tilde. Equivalently, we can normalize a_tilde=1.0

    # Derived parameter: worker separation rate for steady-state
    #s = (τ - δ) / (1 - δ)
end

# =============================================================================
# Equilibrium Functions
# =============================================================================

"""
    H(a)
Productivity-draw cdf given minimal productivity a_m and Pareto shape parameter k 
"""
function H(a, a_m, k)
    (a_m <=0 || k<=0) && throw(ArgumentError("a_m, k must be > 0 "))
    if a >= a_m
        out = 1.0 - (a_m/a)^k 
    else
        out = 0.0 
    end 
    return out
end 

"""
H_inv(cutoff, a_m, k)
Percentile function of productivity draw H. Given cutoff, finds the associated value.
function H_e
"""
function H_inv(cutoff, am, k)
    # Solve cutoff = 1.0 - (a_m/a)^k for a
    (a_m <=0 || k<=0) && throw(ArgumentError("a_m, k must be > 0 "))
    a = a_m/((1-cutoff)^(1/k))
    return a 
end 


"""
    H^e(a)
Survival cdf given cutoff a_star and Pareto shape parameter k 
"""
function H_e(a, a_star, k)
    (a_star <=0 || k<= 0) && throw(ArgumentError("a_star, k must be > 0"))
    if a >= a_star
        out = 1.0 - (a_star/a)^k 
    else
        out = 0.0
    end
    return out 
end

function a_tilde_fun(a_star, k, ε)
    # = E(a^(ε-1)|a>=a*)^(1/(ε-1))
    (k <= (ε-1)) && throw(ArgumentError("Bounds on k violated"))
     Δ = (k/(k-(ε-1)))^(1/(ε-1))
     return Δ*a_star 
end 



"""
    S_p(p, f_r, α)
Top 100p% employment share given ratio of fixed to mean variable labor f_r=f/(l_tilde-f) and shape parameter α 

"""
function S_p(p, f_r, α)
    return p*(f_r + p^(-1/α))/(1+f_r)
end 

function f_r_from_emp_share(α, p::Float64, Sp::Float64)
    Sp = 0.54
    p = 0.02 
    α = 1.1
    f_r = (Sp - p^(1-1/α))/(p-Sp)
    return f_r 
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
    @unpack ρ, δ, F, x_m, ξ_inv = para
    e = e_fun(θ, δ_e, para)
    # Value of vacancy
    Q = (e / F)^ξ_inv * x_m
    return (ρ + δ_e) / (1 + ρ) * Q
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
    θ_fun(para; init_value=0.51)

Solve for steady-state market tightness θ given parameters.
"""
function θ_fun(para; init_value=0.51)
    @unpack f_e, τ, δ, z, b, ϕ, ρ, σ, ε, A, η_L, κ, ξ_inv, F, s = para
    μ = ε / (ε - 1)

    function loss(x)
        θ = abs(x[1]) # ensure positivity
        f = jf(θ, A, η_L)
        q = vf(θ, A, η_L)
        K = K_fun(θ, para)
        L = L_fun(θ, para)
        u = 1 - L
        lhs = (κ + K / q) * (ρ + τ + (1 - δ) * ϕ * q * θ)
        N = (μ - 1) * z * L * (1 - δ) / (f_e * (δ * μ + ρ))
        p = N^(1 / (ε - 1))
        w_int = p * z / μ
        rhs = (1 - δ) * (1 - ϕ) * (w_int - K - b)
        return [(lhs - rhs) / (lhs + rhs)]
    end

    sol = LeastSquaresOptim.optimize(loss, [init_value], Dogleg())
    println("converged=$(sol.converged) at root=$(sol.minimizer) in " *
        "$(sol.iterations) iterations and $(sol.f_calls) function calls")
    return abs(sol.minimizer[1])
end

"""
    steady_state(para; init=0.51)

Compute steady-state statistics given parameters.
Returns a NamedTuple of all key statistics.
"""
function steady_state(para; init=0.51)
    @unpack f_e, τ, δ, z, b, ϕ, ρ, σ, ε, A, η_L, κ, ξ_inv, F, x_m, s = para
    μ = ε / (ε - 1)
    θ = θ_fun(para, init_value=init)

    f = jf(θ, A, η_L)
    q = vf(θ, A, η_L)
    L = L_fun(θ, para)
    u = 1 - L

    # Labor market variables
    u = τ / (τ + (1 - δ) * f)
    v = θ * u
    e = δ * (v + 1 - u)
    Q = (e / F)^ξ_inv * x_m
    K = (ρ + δ) / (1 + ρ) * Q
    X_v = e * (1 / (1 + ξ_inv)) * Q
    v_pret = v - e

    # Relative price, businesses, and values
    N = (μ - 1) * z * L * (1 - δ) / (f_e * (δ * μ + ρ))
    p = N^(1 / (ε - 1))
    N_e = δ / (1 - δ) * N
    ν_f = p * f_e / μ
    d_f = (ρ + δ) / (1 - δ) * ν_f

    # Wages
    w_int = p * z / μ
    w = ϕ * (w_int - K + θ * (K + q * κ)) + (1 - ϕ) * b

    # Sectoral labor 
    L_c = (ρ + δ) * L / (δ * μ + ρ)
    L_e = δ * (μ - 1) * L / (δ * μ + ρ)

    # Consumption output
    Y_c = p * z * L_c

    # Total recruiting costs
    X = X_v + κ * v * q

    # Output and shares
    C = Y_c - X
    Y = C + ν_f * N_e
    J = Q + (1 + ρ) / (1 - δ) * K / q
    M = Q * v + J * L + (N + N_e) * ν_f
    labor_prod = Y / (p * L)

    labor_share = w * L / Y
    cons_share = C / Y
    inv_new_firm_share = ν_f * N_e / Y
    vacancy_share = X / Y
    sunk_vac_cost_share = X_v / Y
    entrant_share = e / v
    x_v = (K / q) / (κ + K / q)
    search_wedge = w / w_int
    recruiter_share = w_int * L / Y
    ann_int_rate = (1 + ρ)^12 - 1

    return (
        θ=θ, N=N, f=f, q=q, u=u, v=v, v_pret=v_pret, e=e, K=K, p=p, N_e=N_e,
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
    x_v=1.0, ξ_inv=1, r_ann=0.04, ε=4.3, σ=1.0, N=1.0, w=1.0, a_m=1.0,
    α=1.1, f_r=0.2)

"""
    calibrate_shares(targets)

Calibrate model parameters to match empirical targets.
Returns a NamedTuple of calibrated parameters.
"""
function calibrate_shares(targets)
    @unpack X_Y, dest_ann, dest_end_frac, f, η_L, q, sep, b_ratio, x_v, ξ_inv, ε, α, r_ann, σ, N, w, a_m = targets

    # Recover Pareto shape of H (productivity draws)
    k = α*(ε-1)

    # Multiplier relating a_tilde and a*
    Δ = (α/(α-1))^(1/(ε-1))

    # Min variable employment = (α-1)/α*(l_tilde-v)

    δ_e = 1 - (1 - dest_ann)^(1 / 12)
    δ = (1.0-dest_end_frac)*δ_e 

    # Purely endogenous destruction threshold
    H_astar = (δ_e-δ)/(1-δ)
    # Back out a_star 
    a_star = H_inv(H_astar, a_m, k)
    a_tilde = Δ*a_star

    @assert H(a_star, a_m, k) ≈ H_astar 

    τ = sep
    s = (τ - δ_e) / (1 - δ_e)
    ρ = (1 + r_ann)^(1 / 12) - 1

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

    #Retailer profit share (would be 1/ε in absence of fixed costs)
    π_s = (1/ε - f_r)



    surplus_ratio = (ρ + τ) / (1 - δ_e) * (1 / (q * x_v))
    p = N^(1 / (ε - 1))
    μ = ε / (ε - 1)
    #recruiter_share = (δ + (ρ + δ) * (ε - 1)) / (δ + (ρ + δ) * ε)
    recruiter_share = ((1-π_s)*ρ+δ_e)/(ρ+δ_e*(1+π_s))

    L_den =  f_e*(δ_e*μ+ρ)+μ*z*f*(1-δ_e)
    L_c = (f_e*(δ_e+ ρ) + μ*z*f*(1-δ_e))*L/L_den
    L_e = f_e*δ_e * (μ - 1) * L / L_den
    @assert L_c + L_e ≈ L 

    function loss(Q)
        Q = abs(Q)
        K = Q * (ρ + δ_e) / (1 + ρ)
        κ = (1 - x_v) / x_v * K / q
        X = e / (1 + ξ_inv) * Q + κ * q * v
        w_int = surplus_ratio * K + w + K
        ϕ = (w - b) / (w_int - K + θ * (K + q * κ) - b)

        z = (μ / p) * w_int *(1/a_tilde) # Now include a_tilde with firm heterogeneity
        #f_e = (μ - 1) * z * (L / N) * (1 - δ) / (δ * μ + ρ)
        ν_f = p * f_e / μ
        Y_c = p * z * L_c
        C = Y_c - X
        Y_new = C + ν_f * N_e
        Y = 1 / X_Y * X
        return (Y - Y_new) / (Y + Y_new), (w_int=w_int, κ=κ, z=z, f_e=f_e, K=K, X=X, Y=Y)
    end

    Q = fzero(x -> loss(x)[1], 1.0)
    out = loss(Q)[2]
    @unpack w_int, κ, z, f_e, K, X, Y = out
    X_Y = abs(X_Y)

    ϕ = (w - b) / (w_int - K + θ * (K + q * κ) - b)
    x_m = Q / e^ξ_inv
    F = 1.0

    return (
        f_e=f_e, τ=τ, δ=δ, z=z, b=b, ϕ=ϕ, ρ=ρ, σ=σ, ε=ε, A=A, η_L=η_L,
        κ=κ, ξ_inv=ξ_inv, x_m=x_m, s=s, F=F
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

