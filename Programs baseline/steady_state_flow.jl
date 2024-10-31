using PyPlot
using Parameters, CSV, StatsBase, Statistics, Random
using NLsolve
using Roots, Optim, LeastSquaresOptim
using PrettyPrinting
using Distributions
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

@with_kw mutable struct ParaCalibFlow

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

    γ::Float64 = 0.4                  # Flow vacancy cost

    # Worker separation rate
    s = (τ - δ)/(1-δ)

end

para = ParaCalibFlow()


function steady_state_flow(para)
    @unpack f_e, τ, δ, z, b, ϕ, ρ, σ, ε, A, η_L, γ, s = para

    μ = ε/(ε-1)

    function loss(x)
        θ = x[1]
        f = jf(θ, A, η_L)
        q = vf(θ, A, η_L)
        L = L_fun(θ, para) 
        u = 1 - L
        lhs = γ*((ρ+τ)/q+ϕ*θ)
        N = (μ-1)*z*L*(1-δ)/(f_e*(δ*μ+ρ))
        p = N^(1/(ε-1))
        w_R = p*z/μ
        rhs = (1-δ)*(1-ϕ)*(w_R-b)
        return [(lhs-rhs)/(lhs+rhs)]
    end

    #θ = fzero(loss, 0.51)
    sol = LeastSquaresOptim.optimize(loss, [0.51], Dogleg())
    println("converged=$(sol.converged) at root=$(sol.minimizer) in " *
        "$(sol.iterations) iterations and $(sol.f_calls) function calls")
    θ = sol.minimizer[1]

    f = jf(θ, A, η_L)
    q = vf(θ, A, η_L)
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
    w_R = p*z/μ
    w = ϕ*(w_R+θ*γ/(1-δ))+(1-ϕ)*(b)

    # Sectoral labor 
    L_e = (δ/(1-δ))*N*f_e/z
    L_c = L-L_e
    # Consumption output
    Y_c = p*z*L_c

    # Additional calculations: V, sunk vacancy costs, C, Y
    C = Y_c - γ*v
    Y = Y_c + ν_f*N_e

    labor_share = w*L/Y
    sunk_entry_cost_share = ν_f*N_e/Y


    out = (θ=θ, N=N, f=f, q=q, u=u, v=v, p=p, N_e=N_e, ν_f=ν_f, d_f=d_f, w_R=w_R, w=w, L=L, L_e=L_e, L_c=L_c, Y_c=Y_c,
     C=C, Y=Y, labor_share=labor_share, sunk_entry_cost_share=sunk_entry_cost_share)
    return out
end


function calibrate_labor_share_flow(targets)
    @unpack labor_share, dest_ann, r_ann, f, η_L, q, sep, b_ratio, ε, σ, N, w = targets

    μ = ε/(ε-1)
    τ = sep
    ρ = (1+r_ann)^(1/12)-1
    β = 1/(1+ρ)
    δ = 1-(1-dest_ann)^(1/12)
    
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

    p = N^(1/(ε-1))
    b = b_ratio*w
    N_e = δ/(1-δ)*N
    recruiter_share = (δ+(ρ+δ)*(ε-1))/(δ+(ρ+δ)*ε) # w_R*L/Y similar to BGM
    w_wR = labor_share/recruiter_share
    w_R = w/(w_wR)
    z = (μ/p)*w_R

    γ = q*(1-δ)/(ρ+τ)*(w_R-w)

    # From wage equation find ϕ
    ϕ = (w-b)/(w_R+θ*γ/(1-δ)-b)

    # N = (μ-1)*zL*(1-δ)/(f_e(δμ+\rho))
    f_e = (μ-1)*z*(L/N)*(1-δ)/(δ*μ+ρ)

    cal = (f_e=f_e, τ=τ, δ=δ, z=z, b=b, ϕ=ϕ, ρ=ρ, σ=σ, ε=ε, A=A, η_L=η_L, γ=γ, s=s)
    return cal
end


function w_fun(θ, N, para)
    @unpack ϕ, b, z, ε, δ, γ = para
    μ = ε/(ε-1)
    w_R = N^(1/(ε-1))*z/μ
    w = (1-ϕ)*b + ϕ*(w_R+θ/(1-δ)*γ)
    return w
end

function L_fun(θ, para)
    @unpack  δ, s, A, η_L = para
    f = jf(θ, A, η_L)
    τ = 1 - (1-δ)*(1-s)
    return (1-δ)*f/(τ+(1-δ)*f)
end

targets = (labor_share=0.66, dest_ann=0.06, r_ann=0.04, f =0.41, η_L=0.6, q=0.8, sep=0.031, b_ratio=0.71, ε=4, σ=1.5, N=1, w=1)
cal = calibrate_labor_share_flow(targets)

steady = steady_state_flow(cal)


@unpack p, L_c, L_e, w, w_R, L, N, N_e, ν_f, d_f, C, Y, Y_c, labor_share, sunk_entry_cost_share = steady
@unpack f_e, τ, δ, z, b, ϕ, ρ, σ, ε, A, η_L, γ, s = cal
μ = ε/(ε-1)
# Accuracy checks
@assert abs(steady.N - targets.N) < 1e-12
@assert abs(steady.f*(1-cal.δ) - targets.f) < 1e-12

@assert abs(N*d_f - Y_c/ε) < 1e-12 # profit share of consumption output
@assert abs(w_R*L/Y - (δ+(ρ+δ)*(ε-1))/(δ+(ρ+δ)*ε)) < 1e-12
@assert abs(Y - w_R*L - N*d_f) < 1e-12
@assert abs(p*z*L_c - w_R*L - N*ν_f*ρ/(1-δ)) < 1e-12
@assert abs(p*z*L_e/μ-ν_f*N_e) < 1e-12
@assert abs(labor_share - (w/w_R)*(δ+(ρ+δ)*(ε-1))/(δ+(ρ+δ)*ε)) < 1e-12

@assert abs(N_e - L_e*z/f_e) < 1e-12

# Steady-state ratios
@show labor_share
@show ϕ
@show sunk_entry_cost_share 
@show C/Y

function N_jcc(θ, para)
    @unpack ρ, τ, A, η_L, δ, ϕ, z, ε, b = para
    q = vf(θ, A, η_L)
    μ = ε/(ε-1)

    w_R = 1/((1-δ)*(1-ϕ))*γ*((ρ+τ)/q+ϕ*θ) + b
    p = w_R*(μ/z)
    N = p^(ε-1)
    return N
end



function N_res(θ, para)
    @unpack  ρ, δ, τ, A, η_L, f_e, ε, z = para
    μ = ε/(ε-1)
    f = jf(θ, A, η_L)
    u = τ/(τ+(1-δ)*f)
    L = 1 - u

    N = (μ-1)*z*L*(1-δ)/(f_e*(δ*μ+ρ))
    return N
end










   