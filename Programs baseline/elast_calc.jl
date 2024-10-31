using PyPlot
using Parameters, CSV, StatsBase, Statistics, Random
using NLsolve
using Roots, Optim, LeastSquaresOptim
using PrettyPrinting
using Distributions
cd(@__DIR__)
include("steady_state.jl")


targets = (labor_share=0.66, dest_ann=0.06, r_ann=0.04, f =0.41, η_L=0.6, q=0.8, sep=0.031, b_ratio=0.71, ξ_inv=1, ε=4, σ=1.5, N=1, w=1)
cal = calibrate_labor_share(targets)
steady = steady_state(cal)

function elast_funcs(para)

    @unpack η_L, ε, δ, τ, ξ_inv, ϕ = para
    μ = ε/(ε-1)
    steady = steady_state(para)
    @unpack θ, u, f, w_R, w, K = steady

    ϵ_e_θ = u*(τ+(1-δ)*f)/(θ*τ+(1-δ)*f)*(θ+(1-η_L)*(1-θ)*(1-u))
    ϵ_K_θ = ξ_inv*ϵ_e_θ

    # argument in term sof ϵ_{θ,z} 
    function loss(x)
        ϵ_wR_z = μ + 1/(ε-1)*(1-η_L)*u*x 
     

        lhs = (ϵ_K_θ+η_L)*x
        rhs = 1/(w_R-w-K)*((1-ϕ)*(w_R*ϵ_wR_z-K*ϵ_K_θ*x)-ϕ*θ*K/(1-δ)*x*(1+ϵ_K_θ))
        return (lhs-rhs)/(lhs+rhs)
    end

    ϵ_θ_z = find_zero(loss, 10.0)
   



