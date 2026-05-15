
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
    @unpack δ, A, η_L = para
    f = jf(θ, A, η_L)
    e = δ*(θ*δ+(1-δ)*f)/(δ+(1-δ)*f)
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
    w = (1-ϕ)*b + ϕ*(z-K+θ*(K+q*κ))
    return w
end

function L_fun(θ, para)
    @unpack  δ, A, η_L = para
    f = jf(θ, A, η_L)
    return (1-δ)*f/(δ+(1-δ)*f)
end

"""
Solve for steady-state θ
"""
function θ_fun(para; init_value=0.51)
    @unpack δ, z, b, ϕ, ρ, A, η_L, ξ_inv, F = para

    function loss(x)
        θ = x[1]
        f = jf(θ, A, η_L)
        q = vf(θ, A, η_L)
        K = K_fun(θ, para)
        L = L_fun(θ, para) 
        u = 1 - L
        lhs = (K/q)*(ρ+δ+(1-δ)ϕ*q*θ)
        # get N from resource constraint curve
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

    # Output
    Y = z*L

    @assert abs(e - δ*(θ*δ+(1-δ)*f)/(δ+(1-δ)*f)) < 1e-12

    # Additional calculations: V, sunk vacancy costs, C, Y

    C = Y - X
    # Stock market cap
    J = Q + (1+ρ)/(1-δ)*K/q
    M = Q*v +J*L 

    labor_share = w*L/Y
    cons_share = C/Y
    vacancy_share = X/Y
    entrant_share = e/v


    out = (θ=θ, f=f, q=q, u=u, v=v, v_pret=v_pret, e=e, K=K, w=w, L=L, Q=Q, X=X, C=C,
     Y=Y, labor_share=labor_share, cons_share=cons_share,
      vacancy_share=vacancy_share, M=M, entrant_share=entrant_share)
    return out
end

function calibrate_labor_share(targets)
    @unpack labor_share, r_ann, f, η_L, q, sep, b_ratio, ξ_inv, w = targets

    δ= sep
    ρ = (1+r_ann)^(1/12)-1
    β = 1/(1+ρ)
    
    # Correct job finding and vacancy filling probablities
    f = f/(1-δ)
    q = q/(1-δ)

    θ = f/q
    u = δ/(δ+(1-δ)*f)
    v = θ*u 

    L = 1 - u
    # Level parameter of matching function 
    A = f/θ^(1-η_L)

    e = δ*(v+1-u)

    b = b_ratio*w
    #w_int = p*z/μ
    z = w/labor_share

    # surplus_ratio = (w_R - w - K)/K 
    surplus_ratio = (ρ+δ)/(1-δ)*(1/(q))
    K = (z-w)/(1+surplus_ratio)

    # From wage equation find ϕ
    ϕ = (w-b)/(z-K+θ*(K)-b)

    # Find F from free entry condition
    #K = (ρ+δ)/(1+ρ)*(e/F)^(1/ξ)
    Q = K*(1+ρ)/(ρ+δ) 
    F = e/Q^(1/ξ_inv)

    cal = (δ=δ, z=z, b=b, ϕ=ϕ, ρ=ρ, A=A, η_L=η_L,κ=κ, ξ_inv=ξ_inv, F=F)
    return cal
end

"""
Calibrate parameters to targets
Targets of N and w reflect choice of units.
N is associated with f_e, and w is associated with z. Normalizing N=1 also implies p=1
"""
function calibrate(targets)
    @unpack ϕ, r_ann, f, η_L, q, sep, b_ratio, ξ_inv, w  = targets

    ρ = (1+r_ann)^(1/12)-1
    β = 1/(1+ρ)
    δ = sep
    
    # Correct job finding and vacancy filling probablities
    f = f/(1-δ)
    q = q/(1-δ)

    θ = f/q
    u = δ/(δ+(1-δ)*f)
    v = θ*u 

    L = 1 - u
    # Level parameter of matching function 
    A = f/θ^(1-η_L)

    e = δ*(v+1-u)

    b = b_ratio*w
    #w_int = p*z/μ

    # surplus_ratio = (w_int - w - K)/K 
    surplus_ratio = (ρ+δ)/(1-δ)*(1/(q))
    K = (1-ϕ)/ϕ*(w-b)/(surplus_ratio +  θ)
    z = surplus_ratio*K + w + K


    # Find F from free entry condition
    #K = (ρ+δ)/(1+ρ)*(e/F)^(1/ξ)
    if ξ_inv > 0
        Q = K*(1+ρ)/(ρ+δ) 
        F = e/Q^(1/ξ_inv)
    end 

    cal = ( δ=δ, z=z, b=b, ϕ=ϕ, ρ=ρ, A=A, η_L=η_L, ξ_inv=ξ_inv, F=F)
    return cal
end

function calibrate_alt(targets)
    """
    Alternative normalization with z = 1 instead of w = 1
    b_ratio is now relative to z, not w, so that b = b_ratio*z = b_ratio
    """

    @unpack ϕ, r_ann, f, η_L, q, sep, b_ratio, ξ_inv, z  = targets

    ρ = (1+r_ann)^(1/12)-1
    β = 1/(1+ρ)
    δ = sep
    
    # Correct job finding and vacancy filling probablities
    f = f/(1-δ)
    q = q/(1-δ)

    θ = f/q
    u = δ/(δ+(1-δ)*f)
    v = θ*u 

    L = 1 - u
    # Level parameter of matching function 
    A = f/θ^(1-η_L)

    e = δ*(v+1-u)

    function loss_fun(w)
        b = b_ratio
        #w_int = p*z/μ

        # surplus_ratio = (w_int - w - K)/K 
        surplus_ratio = (ρ+δ)/(1-δ)*(1/(q))
        K = (1-ϕ)/ϕ*(w-b)/(surplus_ratio +  θ)
        z = surplus_ratio*K + w + K
        loss = z-1
        out = (b=b, K=K, z=z)
        return loss, out
    end

    # Calc wage
    w = fzero(x -> loss_fun(x)[1], 1.0)
    out = loss_fun(w)[2]
    @unpack b, K, z = out

    Q = K*(1+ρ)/(ρ+δ) 
    x_m = Q/e^(ξ_inv)
    #F = e/Q^(1/ξ_inv)
    F = 1

    cal = (δ=δ, z=z, b=b, ϕ=ϕ, ρ=ρ, A=A, η_L=η_L, ξ_inv=ξ_inv, F=F, x_m=x_m)
    return cal
end


function calibration_table(cal, targets)
    @unpack δ, z, b, ϕ, ρ, A, η_L, ξ_inv, A, F = cal
    @unpack ϕ, r_ann, f, η_L, q, sep, b_ratio, ξ_inv, z = targets

    ρ = (1+r_ann)^(1/12)-1
    β = 1/(1+ρ)
    δ = sep

    # Creating a DataFrame for the table
    df = DataFrame(
        Targets = [
            "Real interest rate",
            "Elasticity of matching function",
            "Replacement ratio b/w",
            "Elasticity of vacancy value",
            "Steady-state wage",
            "Aggregate separation rate",
            "Bargaining power",
            "Job finding rate",
            "Vacancy filling rate"
        ],
        Value = round.([r_ann, η_L, b, ξ_inv, w, δ, labor_share, 0.41, 0.80], sigdigits=2),
        Parameter = [L"\rho", L"\eta_L", L"b", L"\delta", L"\xi^{-1}", L"z", L"\phi", L"A", L"F"],
        Calibration = round.([ρ, η_L, b, δ, ξ_inv, z, ϕ, A, F], sigdigits=3)
    )

    # Save the DataFrame as a PDF table
    return df
end






