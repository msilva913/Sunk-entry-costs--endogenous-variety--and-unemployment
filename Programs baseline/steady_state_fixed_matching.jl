
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

    f_e::Float64 = 1                 # Entry cost
    τ::Float64 = 0.0305              # Worker separation rate
    δ::Float64 = 0.005               # Product destruction rate
    z::Float64 = 1.0                    # Technology level
    b::Float64 = 0.71                   # Unemployment insurance
    ϕ::Float64 = 0.5                  # Bargaining power

    ρ::Float64 = 0.04/12                 # Rate of time preference
    σ::Float64 = 1.0                  # Inverse intertemporal elasticity of sub.
    ε::Float64 = 4.0                 # Elasticity of substitution
    
    A::Float64 = 0.5631               # Job matching level parameter
    η_L::Float64 = 0.6                 # Elasticity matching function wrt u

    κ::Float64 = 0.2                 # Fixed matching cost

    # Worker separation rate
    s = (τ - δ)/(1-δ)

end


function w_fun(θ, N, para)
    @unpack ϕ, b, z, ε, δ, A, η_L = para
    q = vf(θ, A, η_L)
    μ = ε/(ε-1)
    w_int = N^(1/(ε-1))*z/μ
    K = K_fun(θ, para)
    w = (1-ϕ)*b + ϕ*(w_int+θ*(q*κ))
    return w
end

function L_fun(θ, para)
    @unpack  δ, s, A, η_L = para
    f = jf(θ, A, η_L)
    τ = 1 - (1-δ)*(1-s)
    return (1-δ)*f/(τ+(1-δ)*f)
end

"""
Solve for steady-state θ
"""
function θ_fun(para; init_value=0.51)
    @unpack f_e, τ, δ, z, b, ϕ, ρ, σ, ε, A, η_L, κ, s = para
    μ = ε/(ε-1)

    function loss(x)
        θ = sqrt(x[1]^2) #ensure positivity
        f = jf(θ, A, η_L)
        q = vf(θ, A, η_L)
        L = L_fun(θ, para) 
        u = 1 - L
        lhs = (κ)*(ρ+τ+(1-δ)ϕ*q*θ)
        # get N from resource constraint curve
        N = (μ-1)*z*L*(1-δ)/(f_e*(δ*μ+ρ))
        p = N^(1/(ε-1))
        w_int = p*z/μ
        rhs = (1-δ)*(1-ϕ)*(w_int-b)
        return [(lhs-rhs)/(lhs+rhs)]
    end

    #θ = fzero(loss, 0.51)
    sol = LeastSquaresOptim.optimize(loss, [init_value], Dogleg())
    println("converged=$(sol.converged) at root=$(sol.minimizer) in " *
        "$(sol.iterations) iterations and $(sol.f_calls) function calls")
    θ = sqrt(sol.minimizer[1]^2)
    return θ
end


function steady_state(para; init=0.51)
    @unpack f_e, τ, δ, z, b, ϕ, ρ, σ, ε, A, η_L, κ, s = para

    μ = ε/(ε-1)
    θ = θ_fun(para)

    f = jf(θ, A, η_L)
    q = vf(θ, A, η_L)
    #K = K_fun(θ, para)
    L = L_fun(θ, para) 
    u = 1 - L

    # Labor market variables
    u = τ/(τ+(1-δ)*f)
    v = θ*u

    # Relative price, businesses, and values
    N = (μ-1)*z*L*(1-δ)/(f_e*(δ*μ+ρ))
    p = N^(1/(ε-1))
    N_e = δ/(1-δ)*N
    ν_f = p*f_e/μ
    d_f = (ρ+δ)/(1-δ)*ν_f

    # Wages
    w_int = p*z/μ
    w = ϕ*(w_int+θ*(q*κ))+(1-ϕ)*(b)

    # Sectoral labor 
    L_c = (ρ+δ)*L/(δ*μ+ρ)
    L_e = δ*(μ-1)*L/(δ*μ+ρ)
    #L_e = (δ/(1-δ))*N*f_e/z
    #L_c = L-L_e
    # Consumption output
    Y_c = p*z*L_c

    # Additional calculations: V, sunk vacancy costs, C, Y
    # total recruiting costs
    X = κ*v*q

    C = Y_c - X
    Y = C + ν_f*N_e
    # Stock market cap
    J = w_int-w + (1-s)*κ
    M = J*L + (N+N_e)*ν_f
    labor_prod = Y/(p*L)

    labor_share = w*L/Y
    cons_share = C/Y
    inv_new_firm_share = ν_f*N_e/Y
    vacancy_share = X/Y
    search_wedge = w/w_int
    recruiter_share = w_int*L/Y
    ann_int_rate = (1+ρ)^12-1

    out = (θ=θ, N=N, f=f, q=q, u=u, v=v, p=p, N_e=N_e, ν_f=ν_f, d_f=d_f, w_int=w_int, w=w,
     L=L, L_e=L_e, L_c=L_c, J=J,
     X=X, C=C, Y_c=Y_c, Y=Y, labor_share=labor_share, 
     labor_prod=labor_prod, cons_share=cons_share, inv_new_firm_share=inv_new_firm_share, vacancy_share=vacancy_share,
      M=M, search_wedge=search_wedge, recruiter_share=recruiter_share, μ=μ, ann_int_rate=ann_int_rate)
    return out
end

targets = (labor_share=0.66, dest_ann=0.0754, f =0.41, η_L=0.6, q=0.8, sep=0.031, b_ratio=0.71, 
         r_ann=0.04, ε=4.3, σ=1.0, N=1.0, w=1.0)


"""
Calibrate parameters to targets
Targets of N and w reflect choice of units.
N is associated with f_e, and w is associated with z. Normalizing N=1 also implies p=1
"""
function calibrate_shares(targets)
    @unpack labor_share, dest_ann, f, η_L, q, sep, b_ratio, ε, r_ann, σ, N, w = targets
    # \nu_f N_e/Y = δ/(ε*(ρ+δ)+δ)
    δ = 1-(1-dest_ann)^(1/12)
    τ = sep
    ρ = (1+r_ann)^(1/12)-1

     # Correct job finding and vacancy filling probablities
    f = f/(1-δ)
    q = q/(1-δ)

    θ = f/q
    u = τ/(τ+(1-δ)*f)
    v = θ*u 

    L = 1 - u
    # Level parameter of matching function 
    A = f/θ^(1-η_L)
    s = (τ-δ)/(1-δ)

    #p = N^(1/(ε-1))
    N_e = δ/(1-δ)*N
    b = b_ratio*w

    # Use labor share and normalization to back out Y 
    Y = w*L/labor_share
    p = N^(1/(ε-1))
    μ = ε/(ε-1) # gross markup 
     # Labor share = (w/w_int)*(recruiter_share)
    recruiter_share = (δ+(ρ+δ)*(ε-1))/(δ+(ρ+δ)*ε) # w_R*L/Y similar to BGM
    L_c = (ρ+δ)*L/(δ*μ+ρ)
    L_e = δ*(μ-1)*L/(δ*μ+ρ)

    # implied recruiting costs 
    function vacancy_loss(X_Y)
        X_Y = sqrt(X_Y^2)
        # Use labor_share = (1+X/Y)*(w/w_int)*recruiter_share
        w_wint = labor_share/(recruiter_share*(1+X_Y))
        w_int = w/(w_wint)
        κ = (1-δ)/(ρ+τ)*(w_int-w)
        X = κ*q*v
        out = (w_int=w_int, κ=κ, X=X)
        return X/Y - X_Y, out
    end

    X_Y = fzero(x -> vacancy_loss(x)[1], 0.015)
    out = vacancy_loss(X_Y)[2]
    @unpack w_int, κ, X = out


    #w_int = p*z/μ
    z = (μ/p)*w_int
    # Consumption output
    Y_c = p*z*L_c

    # From wage equation find ϕ
    ϕ = (w-b)/(w_int+θ*(q*κ)-b)
    # Sectoral labor  
    f_e = (μ-1)*z*(L/N)*(1-δ)/(δ*μ+ρ)
    ν_f = p*f_e/μ

    cal = (f_e=f_e, τ=τ, δ=δ, z=z, b=b, ϕ=ϕ, ρ=ρ, σ=σ, ε=ε, A=A, η_L=η_L,κ=κ, s=s)

    return cal
end


"""
Job creation curve:
N as a function of θ
"""
function N_jcc(θ, para)
    @unpack ρ, τ, A, η_L, δ, ϕ, z, ε, b, κ = para
    q = vf(θ, A, η_L)
    μ = ε/(ε-1)

    w_int = 1/((1-δ)*(1-ϕ))*(κ*q)*((ρ+τ)/q+(1-δ)*ϕ*θ) + b
    p = w_int*(μ/z)
    N = p^(ε-1)
    return N
end


"""
Resource constraint curve:
N as a function of θ
"""
function N_res(θ, para)
    @unpack  ρ, δ, τ, A, η_L, f_e, ε, z = para
    μ = ε/(ε-1)
    f = jf(θ, A, η_L)
    u = τ/(τ+(1-δ)*f)
    L = 1 - u

    N = (μ-1)*z*L*(1-δ)/(f_e*(δ*μ+ρ))
    return N
end






