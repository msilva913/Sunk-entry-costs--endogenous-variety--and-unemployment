
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

    F::Float64 = 0.0007              # Mass of recruiters with opportunity to create vacancy
    κ::Float64 = 0.2                 # Fixed matching cost
    ξ_inv::Float64 = 1.0             # Inverse elasticity of entry to vacancy value

    # Worker separation rate
    s = (τ - δ)/(1-δ)

end

# Equilibrium curves
function e_fun(θ, para)
    # e as a function of θ
    @unpack δ, τ, A, η_L = para
    f = jf(θ, A, η_L)
    e = δ*(θ*τ+(1-δ)*f)/(τ+(1-δ)*f)
    return e
end

function K_fun(θ, para)
    @unpack ρ, δ, F, ξ_inv = para
    e = e_fun(θ, para)
    K = (ρ+δ)/(1+ρ)*(e/F)^(ξ_inv)
    return K
end

function w_fun(θ, N, para)
    @unpack ϕ, b, z, ε, δ, A, η_L = para
    q = vf(θ, A, η_L)
    μ = ε/(ε-1)
    w_int = N^(1/(ε-1))*z/μ
    K = K_fun(θ, para)
    w = (1-ϕ)*b + ϕ*(w_int-K+θ*(K+q*κ))
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
    @unpack f_e, τ, δ, z, b, ϕ, ρ, σ, ε, A, η_L, κ, ξ_inv, F, s = para
    μ = ε/(ε-1)

    function loss(x)
        θ = x[1]
        f = jf(θ, A, η_L)
        q = vf(θ, A, η_L)
        K = K_fun(θ, para)
        L = L_fun(θ, para) 
        u = 1 - L
        lhs = (κ+K/q)*(ρ+τ+(1-δ)ϕ*q*θ)
        # get N from resource constraint curve
        N = (μ-1)*z*L*(1-δ)/(f_e*(δ*μ+ρ))
        p = N^(1/(ε-1))
        w_int = p*z/μ
        rhs = (1-δ)*(1-ϕ)*(w_int-K-b)
        return [(lhs-rhs)/(lhs+rhs)]
    end

    #θ = fzero(loss, 0.51)
    sol = LeastSquaresOptim.optimize(loss, [init_value], Dogleg())
    println("converged=$(sol.converged) at root=$(sol.minimizer) in " *
        "$(sol.iterations) iterations and $(sol.f_calls) function calls")
    θ = sol.minimizer[1]
    return θ
end



function calibrate_labor_share(targets)
    @unpack labor_share, dest_ann, r_ann, f, η_L, q, sep, b_ratio, x_v, ξ_inv, ε, σ, N, w = targets

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

    e = δ*(v+1-u)

    p = N^(1/(ε-1))
    N_e = δ/(1-δ)*N
    b = b_ratio*w
    #w_int = p*z/μ
    recruiter_share = (δ+(ρ+δ)*(ε-1))/(δ+(ρ+δ)*ε) # w_R*L/Y similar to BGM
    w_wint = labor_share/recruiter_share
    w_int = w/(w_wint)
    z = (μ/p)*w_int

    # surplus_ratio = (w_R - w - K)/K 
    surplus_ratio = (ρ+τ)/(1-δ)*(1/(q*x_v))
    K = (w_int-w)/(1+surplus_ratio)

    # Find κ given K 
    κ = (1-x_v)/x_v*K/q

    # From wage equation find ϕ
    ϕ = (w-b)/(w_int-K+θ*(K+q*κ)-b)

    # Find F from free entry condition
    #K = (ρ+δ)/(1+ρ)*(e/F)^(1/ξ)
 
    Q = K*(1+ρ)/(ρ+δ) 
    F = e/Q^(1/ξ_inv)
    # Given N, solve for f_e
    # N = (μ-1)*zL*(1-δ)/(f_e(δμ+\rho))
    f_e = (μ-1)*z*(L/N)*(1-δ)/(δ*μ+ρ)

    cal = (f_e=f_e, τ=τ, δ=δ, z=z, b=b, ϕ=ϕ, ρ=ρ, σ=σ, ε=ε, A=A, η_L=η_L,κ=κ, ξ_inv=ξ_inv, F=F, K=K, s=s)
    return cal
end

Traceplot_parameter = matopen("Traceplot_parameter.mat")
Traceplot_parameter = read(Traceplot_parameter, "Traceplot_parameter")
# For steady-state-based dependent parameters, we can omit Shocks
Traceplot_parameter = Traceplot_parameter[:, 1:(end-6)]
col_names = [:σ, :b, :x_v, :ξ_inv, :δ, :ε]
df = DataFrame(Traceplot_parameter, :auto)
rename!(df, col_names)

#labor_share, dest_ann, r_ann, f, η_L, q, sep, b_ratio, x_v, ξ_inv, ε, σ, N, w = targets
#targets = (labor_share=0.66, dest_ann=0.10, r_ann=0.04, f =0.41, η_L=0.6, q=0.8, sep=0.031, b_ratio=0.71, x_v=0.1, ξ_inv=1/0.265, ε=4.3, σ=1.0, N=1.0, w=1.0)
n = size(df, 1)
targets_df = DataFrame(
labor_share = 0.66*ones(n),
dest_ann = 1.0 .-((1. .-df[!,:δ]).^12),
r_ann = 0.04*ones(n),
f = 0.41*ones(n),
η_L = 0.6*ones(n),
q = 0.8*ones(n),
sep = 0.031*ones(n), 
b_ratio = df[!, :b],
x_v = df[!, :x_v],
ξ_inv = df[!, :ξ_inv],
ε = df[!, :ε],
σ = df[!, :σ],
N = 1.0*ones(n),
w = 1.0*ones(n),
)

results = [calibrate_labor_share(targets_df[row,:]) for row in 1:n]

dep_parameters_df = DataFrame(results)

κ = dep_parameters_df[!, :κ]
ϕ = dep_parameters_df[!,:ϕ]
K = dep_parameters_df[!,:K]
using KernelDensity
filtered_data(data) = filter(!isnan, data)

density = kde(filtered_data(ϕ))

fig, ax = plt.subplots()
ax.plot(density.x, density.density, color="orange", alpha=0.3)
ax.fill_between(density.x, density.density, color="gold", alpha=0.3, zorder=1)
ax.set_xlabel("ϕ", fontsize=12)
ax.set_ylabel("Density", fontsize=12)
display(fig)

fig, ax = plt.subplots()
ax.hist(K)
display(fig)

cd("C:/Users/msilva913/Documents/GitHub/Sunk_entry_costs_endogenous_variety_unemployment/Programs baseline/")
include(pwd()*"/steady_state.jl")
include(pwd()*"/generate_tables.jl")

cd(@__DIR__)
posterior_mean= matopen("posterior_mean.mat")
posterior_mean = read(posterior_mean, "posterior_mean")

delta_mean = posterior_mean["delta"]
delta_ann_mean = 1 - (1-delta_mean)^12

targets = (labor_share=0.66, 
           dest_ann=delta_ann_mean, 
           r_ann=0.04, 
           f =0.41, 
           η_L=0.6, #elast. of matching function
           q=0.8, 
           sep=0.031, 
           b_ratio=posterior_mean["b_ratio"], 
            x_v=posterior_mean["xi_inv"], 
            ξ_inv=posterior_mean["xi_inv"], # congestion elasticity
            ε=posterior_mean["epsi"], # elasticity of sub.
            σ=posterior_mean["sigma"], # log utility
            N=1, w=1.0)

cal = calibrate_labor_share(targets)
ss = steady_state(cal)
cal_table = calibration_table(cal, targets)

output = IOBuffer()
show(output, MIME("text/latex"),cal_table)
cal_table_tex = String(take!(output))
print(cal_table_tex)

df_shares = shares_table(ss)

output = IOBuffer()
show(output, MIME("text/latex"),df_shares)
df_shares_table = String(take!(output))
print(df_shares_table)
