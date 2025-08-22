
using Parameters, CSV, StatsBase, Statistics, Random
using NLsolve
using DataFrames
using Roots, Optim, LeastSquaresOptim
using PrettyPrinting
using LaTeXStrings
cd(@__DIR__)
# Matching probabilities
function jf(θ, A, η_L)
    return min(A*θ^(1-η_L), 1.0)
end

#vacancy filling probability
function vf(θ, A, η_L)
    return min(A*θ^(-η_L), 1.0)
end

function θ_invert(q, A, η_L)
    θ = (q/A)^(-1/η_L)
    return θ
end

@with_kw mutable struct ParaCalib

    δ::Float64 = 0.005               # Firm destruction rate/worker separation rate
    z::Float64 = 1.0                    # Technology level
    b::Float64 = 0.71                   # Unemployment insurance
    ϕ::Float64 = 0.5                  # Bargaining power

    ρ::Float64 = 0.04/12                 # Rate of time preference
   
    A::Float64 = 0.5631               # Job matching level parameter
    η_L::Float64 = 0.6                 # Elasticity matching function wrt u

    F::Float64 = 0.0007              # Mass of recruiters with opportunity to create vacancy
    ξ_inv::Float64 = 1.0             # Inverse elasticity of entry to vacancy value
    x_m::Float64 = 10                # Level parameter of cost dist.


end

# Equilibrium curves
function e_fun(θ, para)
    # e as a function of θ
    @unpack δ, s, A, η_L = para
    τ=s+δ*(1-s)
    f = jf(θ, A, η_L)
    e = δ * (θ * τ + (1 - δ) * f) / (τ + (1 - δ) * f)
    return e
end

function K_fun(θ, para)
    @unpack ρ, δ, F, x_m, ξ_inv = para
    e = e_fun(θ, para)
    Q = (e/F)^ξ_inv*x_m
    K = (ρ+δ)/(1+ρ)*Q
    return K
end

function w_fun(θ, N, para)
    @unpack ϕ, b, z, δ, A, η_L = para
    q = vf(θ, A, η_L)
    K = K_fun(θ, para)
    w = (1-ϕ)*b + ϕ*(z-K+θ*(K))
    return w
end

function L_fun(θ, para)
    @unpack  δ, s, A, η_L = para
    τ = s + δ*(1-s)
    f = jf(θ, A, η_L)
    return (1-δ)*f/(τ+(1-δ)*f)
end

"""
Solve for steady-state θ
"""
function θ_fun(para; init_value=0.51)
    @unpack δ, s, z, b, ϕ, ρ, A, η_L, ξ_inv, F, x_m = para

    τ = s + δ*(1-s)

    function loss(x)
        θ = x[1]
        f = jf(θ, A, η_L)
        q = vf(θ, A, η_L)
        K = K_fun(θ, para)
        L = L_fun(θ, para) 
        u = 1 - L
        lhs = (K/q)*(ρ+ τ +(1-δ)*ϕ*q*θ)
        rhs = (1-δ)*(1-ϕ)*(z-K-b)
        return [(lhs-rhs)/(lhs+rhs)]
    end

    #θ = fzero(loss, 0.51)
    sol = LeastSquaresOptim.optimize(loss, [init_value], Dogleg())
    println("converged=$(sol.converged) at root=$(sol.minimizer) in " *
        "$(sol.iterations) iterations and $(sol.f_calls) function calls")
    θ = sol.minimizer[1]
    return θ
end


function steady_state(para; init=0.51)
    @unpack δ, z, b, ϕ, ρ, A, η_L, ξ_inv, F, x_m = para

    θ = θ_fun(para)

    f = jf(θ, A, η_L)
    q = vf(θ, A, η_L)
    K = K_fun(θ, para)
    L = L_fun(θ, para) 
    u = 1 - L

    # Labor market variables
    u = δ/(δ+(1-δ)*f)
    v = θ*u
    e = δ*(v+1-u)
    Q = (e/F)^ξ_inv*x_m
    K = (ρ+δ)/(1+ρ)*Q
    X = e/(1+ξ_inv)*Q
   
    v_pret = v-e

    w = ϕ*(z-K+θ*K)+(1-ϕ)*(b)

    # Gross Output
    Y = z*L

    @assert abs(e - δ*(θ*δ+(1-δ)*f)/(δ+(1-δ)*f)) < 1e-12

    # Additional calculations: V, sunk vacancy costs, C, Y

    C = Y - X
    # Stock market cap
    J = Q + (1+ρ)/(1-δ)*K/q
    M = Q*v +J*L 

    vacancy_cost_share = X/C
    entrant_share = e/v
    profit_share = ((z-w)*(1-u) - X)/C
    labor_share = w*(1-u)/C


    out = (θ=θ, f=f, q=q, u=u, v=v, v_pret=v_pret, e=e, K=K, w=w, L=L, Q=Q, X=X, C=C,
     Y=Y, vacancy_cost_share=vacancy_cost_share, M=M, entrant_share=entrant_share, profit_share=profit_share,
     labor_share=labor_share)
    return out
end

targets = (
    X_Y=0.015, dest_ann=0.0754, f=0.41, η_L=0.6, q=0.8, sep=0.031, b_ratio=0.71,
    ξ_inv=1, r_ann=0.04, w=1.0
)

function calibrate_shares(targets)
    @unpack X_Y, dest_ann, f, η_L, q, sep, b_ratio, ξ_inv, r_ann, w = targets

    # Separation rates
    δ = 1 - (1 - dest_ann)^(1 / 12)
    τ = sep
    s = (τ - δ) / (1 - δ)

    # Rate of time preference
    ρ = (1 + r_ann)^(1 / 12) - 1
    
    # Correct job finding and vacancy filling probablities
    f = f/(1-δ)
    q = q/(1-δ)

    θ = f/q
    u = τ/(τ+(1-δ)*f)
    v = θ*u 

    # Level parameter of matching function 
    A = f/θ^(1-η_L)
    e = δ * (v + 1 - u)

    b = b_ratio*w
    
    # surplus_ratio = (z - w - K)/K 
    surplus_ratio = (ρ+ τ)/(1-δ)*(1/(q))

     function loss(Q)
        Q = abs(Q)
        K = Q * (ρ + δ) / (1 + ρ)
        X = e / (1 + ξ_inv) * Q 
        z = surplus_ratio * K + w + K
        Y_new = z * (1-u) - X
        Y = 1 / X_Y * X
        return (Y - Y_new) / (Y + Y_new), (z=z, K=K, X=X, Y=Y)
    end

    Q = fzero(x -> loss(x)[1], 1.0)
    out = loss(Q)[2]
    @unpack z, K, X, Y = out

    # From wage equation find ϕ
    ϕ = (w-b)/(z-K+θ*(K)-b)
    x_m = Q / e^ξ_inv
    F = 1.0

    cal = (δ=δ, s=s, z=z, b=b, ϕ=ϕ, ρ=ρ, A=A, η_L=η_L, ξ_inv=ξ_inv, x_m=x_m, F=F)
    return cal
end

targets = (
    X_Y=0.015, dest_ann=0.0754, f=0.41, η_L=0.6, q=0.8, sep=0.031, b_ratio=0.71,
    ξ_inv=1, r_ann=0.04, w=1.0
)

function calibration_table(cal, targets)
    @unpack δ, s, z, b, ϕ, ρ, A, η_L, ξ_inv, A, x_m, F = cal
    @unpack X_Y, dest_ann, f, q, sep, b_ratio,  ξ_inv, r_ann, w  = targets


    β = 1/(1+ρ)

    # Creating a DataFrame for the table
    df = DataFrame(
        Targets = [
            "Real interest rate",
            "Product destruction rate",
            "Elasticity of matching function",
            "Replacement ratio b/w",
            "Elasticity of vacancy value",
            "Steady-state wage",
            "Aggregate separation rate",
            "Recruiting cost share",
            "Job finding rate",
            "Vacancy filling rate",
            "Mass of recruiters"
        ],
        Value = round.([r_ann, δ, η_L, b, ξ_inv, w, sep, X_Y, f, q, F], sigdigits=2),
        Parameter = [L"\rho", L"\delta", L"\eta_L", L"b", L"\xi^{-1}", L"z", L"s", L"\phi", L"A", L"x_m", L"F"],
        Calibration = round.([ρ, δ, η_L, b, ξ_inv, z, s, ϕ, A, x_m, F], sigdigits=3)
    )

    # Save the DataFrame as a PDF table
    return df
end

#ss = steady_state(cal)





