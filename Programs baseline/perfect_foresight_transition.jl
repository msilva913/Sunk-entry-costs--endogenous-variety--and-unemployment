using PyPlot
using Parameters, CSV, StatsBase, Statistics, Random
using Roots, Optim, LeastSquaresOptim, NLsolve
using Printf
using PrettyPrinting
using DataFrames
cd(@__DIR__)
include("steady_state.jl")


targets = (labor_share=0.66, dest_ann=0.06, r_ann=0.04, f =0.41, η_L=0.6, q=0.8, sep=0.031, b_ratio=0.71, ξ=1, ε=4, σ=1.5, N=1, z=1)
cal = calibrate_labor_share(targets)
para = cal
steady = steady_state(cal)

@with_kw mutable struct Model_Arrays
    T::Int64 = 50
    v::Vector{Float64} = zeros(T)
    N::Vector{Float64} = zeros(T)
    N_e::Vector{Float64} = zeros(T)
    u::Vector{Float64} = zeros(T)
    θ::Vector{Float64} = zeros(T)
    f::Vector{Float64} = zeros(T)
    q::Vector{Float64} = zeros(T)
    v_pret::Vector{Float64} = zeros(T)
    e::Vector{Float64} = zeros(T)
    p::Vector{Float64} = zeros(T)
    w_R::Vector{Float64} = zeros(T)
    w::Vector{Float64} = zeros(T)
    ν_f::Vector{Float64} = zeros(T)
    Q::Vector{Float64} = zeros(T)
    Ω::Vector{Float64} = zeros(T)
    K::Vector{Float64} = zeros(T)
    L_c::Vector{Float64} = zeros(T)
    L_e::Vector{Float64} = zeros(T)
    C::Vector{Float64} = zeros(T)
    λ::Vector{Float64} = zeros(T)
    Y_c::Vector{Float64} = zeros(T)
    d_f::Vector{Float64} = zeros(T)
    K_new::Vector{Float64} = zeros(T)
    ν_f_new::Vector{Float64} = zeros(T)
end

m = Model_Arrays(T=50)
@unpack v, N_e = m 
v .= steady.v
N_e .= steady.N_e
m.v, m.N_e = v, N_e

function construct_sequences(para, m, shock_seq)
    @unpack  f_e, τ, δ, z, b, ϕ, ρ, σ, ε, A, η_L, ξ, F, s = para
    @unpack N, v, p, w_R, ν_f, d_f, N_e, L_e, L_c, Y_c, T = m
    @unpack θ, f, q, u, v_pret, e, Q, Ω = m
    @unpack Y_c, Ω, C, λ, w, w_R, T = m 

    z_seq = shock_seq

    # Construct number of firms 
    for t in 2:T
        N[t] = (1-δ)*(N[t-1]+N_e[t-1]) #N_2 is predetermined since it depends only on N_1 and N_e1
    end

    μ = ε/(ε-1.0)
    β = 1/(1+ρ)

    @. p .= N^(1/ε)
    @. w_R .= p*z/μ
    @. ν_f .= p*f_e/μ



    for t in 1:(T)
        t1 = min(t+1, T)
        θ[t] = v[t]/u[t]
        f[t] = jf(θ[t], A, η_L)
        q[t] = vf(θ[t], A, η_L)
        u[t1] = (1-(1-δ)*f[t])*u[t] + τ*(1-u[t])
    end

    # Calculate predetermined vacancies and entrants
    for t in 2:T
        v_pret[t] = (1-δ)*((1-q[t-1])*v[t-1]+s*(1-u[t-1]))
        v_pret[t] = min(v_pret[t], v[t])
        e[t] = v[t] - v_pret[t]
        Q[t] = (e[t]/F)^(1/ξ)
        Ω[t] = F*ξ/(ξ+1)*Q[t]^(ξ+1)
    end

    # Labor
    L = @.  1.0- u
    @. L_e .= N_e*f_e/z_seq
    @. L_c .= L - L_e

    # Retail output Y_c and profits 
    @. Y_c .= p*z_seq*L_c
    @. d_f .= Y_c/(N*ε)

    C .= Y_c - Ω
    C .= max.(C, 1e-6)
    @. λ .= C^(-σ)

    for t = 2:T
        t1 = min(t+1,T)
        K[t] = Q[t] - β*(1-δ)*λ[t1]/λ[t]*Q[t1]
    end


    @. w = ϕ*(w_R -K+θ/(1-δ)*K) + (1-ϕ)*b

    m.u, m.θ, m.f, m.q = u, θ, f, q
    m.p, m.w_R, m.ν_f, m.N, m.L_e, m.L_c, m.Y_c, m.d_f = p, w_R, ν_f, N, L_e, L_c, Y_c, d_f
    m.v_pret, m.e, m.Q, m.Ω = v_pret, e, Q, Ω
    m.C, m.λ, m.K, m.w = C, λ, K, w 
    return m
end

function get_updates(para, m, shock_seq, damp=0.5)
    @unpack  f_e, τ, δ, z, b, ϕ, ρ, σ, ε, A, η_L, ξ, F, s = para
    @unpack q, λ, w_R, w, ν_f, d_f, Q, e, v_pret, L_e, K_new, ν_f_new, T = m

    z_seq = shock_seq

    β = 1/(1+ρ)
    μ = ε/(ε-1)

    # Start at T and work backwards to 2
    for t in (T):-1:2
        t1 = min(t+1, T)
        K_new[t] = q[t]*β*(λ[t1]/λ[t])*(1-δ)*(w_R[t1]-w[t1]-K[t1]+(1-s)*K[t1]/q[t1])
        ν_f_new[t] =β*(1-δ)*(λ[t1]/λ[t])*(ν_f[t1]+d_f[t1])
    end

    ii = 2:T 
    loss = zeros(T-1, 2)
    loss .=  hcat((K_new[ii] - K[ii])./K[ii], (ν_f_new[ii]-ν_f[ii])./ν_f[ii])

    @. K_new .= damp*K_new + (1-damp)*K
    @. ν_f_new .= damp*ν_f_new + (1-damp)*ν_f

    # Update N_e

    p = @. ν_f_new*μ/f_e
    N = @. p^(ε-1)

    for t in 2:T
        t1 = min(t+1, T)
        N_e[t] =N[t1]/(1-δ) - N[t]
        N_e[t] = max(N_e[t], 0.0)
        # Enforce bound on entrants from labor supply
        N_e[t] = min(N_e[t], L_e[t]*z_seq[t]/f_e)
    end



    # Update v 

    # Back out Q
    for t in T:-1:2
        t1 = min(t+1, T)
        Q[t] = K_new[t] + β*(1-δ)*(λ[t1]/λ[t])Q[t1]
    end

    # Obtain e from Q 
    e = @. F*Q^(ξ)
    v = @. v_pret + e
    m.N_e = N_e
    m.v = v
    return m, loss
end

function cons_shock_seq(para, m; scal=1.0, ρ_z=0.979)
    @unpack z, δ = para 
    @unpack T = m
    z_seq = zeros(T)
    z_seq[1] = 0
    z_seq[2] =scal/100.0
    # Shock in log deviations
    for t in 3:T
        z_seq[t] = ρ_z*z_seq[t-1]
    end
    # Exponentiate and scale by parameter value
    z_seq .= exp.(z_seq)*z

    shock_seq = (z_seq=z_seq)
    return shock_seq
end

function get_transition!(para, m; tol=1e-5, itermax=1000, T=100)
   
    # Set arrays to steady state initially
    m = Model_Arrays(T=T)
    @unpack v, N, N_e, u, θ, f, q, v_pret, e, p, w_R, w, ν_f, Q, Ω, K, L_c, L_e, C, Y_c, d_f, K_new, ν_f_new = m
    v .= steady.v
    N_e .= steady.N_e
    m.v, m.N_e = v, N_e

    N[1] = steady.N
    u[1] = steady.u 
    θ[1] = steady.θ
    f[1] = steady.f
    q[1] = steady.q 
    v_pret[1] = steady.v_pret[1] 
    e[1] = steady.e 
    p[1] = steady.p 
    w_R[1] = steady.w_R 
    w[1] = steady.w 
    ν_f[1] = steady.ν_f 
    Q[1] = steady.Q 
    Ω[1] = steady.Ω 
    K[1] = steady.K 
    L_c[1] = steady.L_c 
    L_e[1] = steady.L_e 
    C[1] = steady.C 
    Y_c[1] = steady.Y_c 
    d_f[1] = steady.d_f 
    K_new[1] = steady.K 
    ν_f_new[1] = steady.ν_f

    m.N, m.u, m.θ, m.f, m.q, m.v_pret, m.e, m.p, m.w_R, m.w, m.ν_f, m.Q, m.Ω, m.K, m.L_c, m.L_e, m.C, m.Y_c, m.d_f, m.K_new, m.ν_f_new =
    (N, u, θ, f, q, v_pret, e, p, w_R, w, ν_f, Q, Ω, K, L_c, L_e, C, Y_c, d_f, K_new, ν_f_new)

    shock = cons_shock_seq(para, m, scal=1.0, ρ_z=0.979)

    iter = 1
    max_err = 1.0
    while max_err > tol && (iter <= itermax)
    #while iter <= itermax
        m = construct_sequences(para, m, shock)
        m, loss = get_updates(para, m, shock, 0.2)

        # check for number of markets in equilibrium
        max_err = maximum(abs.(loss))
        @printf "Max error in market clearing %.04f\n" max_err
        iter += 1
    end

    if (iter <= itermax)
        @printf "Converged in %iterations " iter
    else
        @printf "No equilibrium found"
    end
end










function factor_prices!(it, para, m)
    # calculate factor prices in a given year
    @unpack γ, egam, β, α, δ, κ, g = para
    @unpack K, L, r, τ_w, τ_p, τ_c, τ_r, w, wn, Rn, p, pen = m
    r[it] = α*(K[it]/L[it])^(α-1) - δ
    w[it] = (1-α)*(K[it]/L[it])^(α)
    wn[it] = w[it]*(1-τ_w[it]-τ_p[it])
    Rn[it] = 1.0 +r[it]*(1-τ_r[it])
    p[it] = 1.0 + τ_c[it]
    pen[it] = κ*w[max(it-1, 1)]
    m.r, m.w, m.wn, m.Rn, m.p, m.pen = r, w, wn, Rn, p, pen
    return m
end

function year(it, ij, ijp)
    year = it + (ijp - ij)
    if (it == 1 || year <= 1)
        year = 1
    elseif (it == T || year >= T)
        year = T 
    end
    return year
end


function decisions!(it, para)
    @unpack γ, egam, β, α, δ, κ, g = para
    # calculate years
    it1 = year(it, 1, 2)
    it2 = year(it, 1, 3)
    itm = year(it, 2, 1)


    PVI = wn[it] + wn[it1]/Rn[it1] + pen[it2]/(Rn[it1]*Rn[it2]) + v[it]
    Ψ = (1/p[it])*(1.0 + β^γ*(p[it1]/(p[it]*Rn[it1]))^(1-γ) + β^(2γ)*(p[it2]/(p[it]*Rn[it1]*Rn[it2]))^(1-γ))^(-1)
    
    c[1, it] = PVI*Ψ
    #it =2: first transition period--first period after initial steady state
    if it == 2
        PVI = Rn[it]*a[2, 1] + wn[it] + pen[it1]/Rn[it1] + v[2]
        Ψ = (1/p[it])*(1.0 + β^γ*(p[it1]/(p[it]*Rn[it1]))^(1-γ))^(-1)
        c[2, it] = PVI*Ψ
        c[3, it] = (pen[it] + Rn[it]*a[3, itm] + v[1])/p[it]
        a[2, it] = wn[itm] - p[itm]*c[1, itm]
    else
        c[2, it] = (β*Rn[it]*p[itm]/p[it])^γ*c[1, itm]
        c[3, it] = (β*Rn[it]*p[itm]/p[it])^γ*c[2, itm]
        a[2, it] = wn[itm] - p[itm]*c[1, itm]
    end
    # initial middle aged who are now old
    if it == 3
        a[3, it] = wn[itm] +  Rn[itm]*a[2, itm] + v[2] - p[itm]*c[2, itm]
    else
        a[3, it] = wn[itm] + Rn[itm]*a[2, itm] - p[itm]*c[2, itm]
    end
    return nothing
end

function utility_fun!(it, para)
    @unpack β, egam = para
    # get future years
    # at steady state
    if it == 1
        it1 = 1
        it2 = 1
    else
        it1 = min(it+1, T)
        it2 = min(it+2, T)
    end

    # oldest cohort
    util[3, it] = c[3, it]^(egam)/egam
    # middle cohort
    util[2, it] = c[2, it]^(egam)/egam + β*c[3, it1]^egam/egam
    # youngest cohort
    util[1, it] = c[1, it]^egam/egam + β*c[2, it1]^(egam)/egam + β^2*c[3, it2]^(egam)/egam
    return nothing
end


# function lsra!(PVI, BA, v, para)
#     @unpack γ, egam, β, α, δ, by, κ, n_p, g = para
#     PVI = zeros(T+1)
#     util = zeros(3, T)
#     for it in 2:T 
#         utility_fun(util, c, it, para)
#     end
#     # Transfers to old generations
#     PVI[1] = Rn[2]*a[3, 2] + pen[2] + v[1]
#     PVI[2] = Rn[2]*a[2, 2] + wn[2] + pen[3]/Rn[3] + v[2]
#     v[1] = v[1] + PV[1]*((util[3, 1]/util[3, 2])^(1/egam) - 1)
#     v[2] = v[2] + PV[2]*((util[2, 1]/util[2, 2])^(1/egam) - 1)
#     BA[3] = v[1]/(1+n_p[1])*(1+n_p[2]) + v[1]/(1+n_p[2])

#     # long-run equilibrium
#     PVI[T+1] = wn[T] + wn[T]/Rn[T] + pen[T]/Rn[T]^2 + v[T+1]
#     PV = v[T+1]*(1.0+r[T])/(r[T]-n_p[T])
#     sum1 = PVI[T+1]*(1.0+r[T])/(r[T]-n_p[T])
#     sum2 = PVI[T+1]*(util[2, TT]*egam)^(-1/egam)*(1+r[T])/(r[T]-n_p[T])

#     # transition path
#     for it in (T-1):-1:2
#         it1 = year(it, 1, 2)
#         it2 = year(it, 1, 3)
#         PVI[it+1] = wn[it] + wn[it1]/Rn[it1] + pen[it2]/(Rn[it1]*Rn[it2]) + v[it+1]
#         PV = PV*(1+n_p[it1])/(1+r[it1]) + v[it+1]
#         sum1 = sum1*(1+n_p[it1])/(1+r[it1]) + PVI[it+1]
#         sum2 = sum2*(1+n_p[it1])/(1+r[it1]) + PVI[it+1]*(util[1, it]*egam)^(-1/egam)
#     end

#     # calculate ustar for future generations
#     ustar = ((sum1-BA[2]-PV)/sum2)^(egam)/egam

#     # calculate transfers to future generations and ebt of LSRA
#     for it in 2:T 
#         v[it+1] = v[it+1] + PVI[it+1]*((ustar/util[1, it])^(1/egam) - 1.0)
#         if it == 3
#             BA[3] = (BA[3] + v[2])/(1+n_p[2])
#         elseif it > 3
#             BA[it] = (1.0+r[it-1])*BA[it-1] + v[it]/(1+n_p[it])
#     end
# end

function quantities!(it, para; damp=0.25)
    @unpack γ, egam, β, α, κ, g, δ = para
    itm = year(it, 2, 1)
    Y[it] = K[it]^α*L[it]^(1-α)
    C[it] = c[1, it] + c[2, it]/(1+n_p[it]) + c[3, it]/((1+n_p[it])*(1+n_p[itm]))
    G[it] = g[1] + g[2]/(1+n_p[it]) + g[3]/((1+n_p[it])*(1+n_p[itm]))
    A[it] = a[2, it]/(1+n_p[it]) + a[3, it]/((1+n_p[it])*(1+n_p[itm]))
    B[it] = by[itm]*Y[it]

    # starting from t=2: after steady state
    if it > 1
        # A = K + B
        K[it] = damp*(A[it]-B[it]-BA[it]) + (1-damp)*K[it]
        it1 = min(it+1, T)
        I[it] = (1+n_p[it1])*K[it1] - (1-δ)*K[it]
    end
end

function government!(it, para)
    it1 = min(it+1, T)
    if tax[it] == 1
        τ_c[it] = ((1+r[it])*B[it] + G[it] - (τ_w[it]*w[it]*L[it] + τ_r[it]*r[it]*A[it] +
                     (1+n_p[it])*B[it1]))/C[it]
    elseif tax[it] == 2
        τ_w[it] = ((1+r[it])*B[it] + G[it] - (τ_c[it]*C[it] + (1+n_p[it])*B[it1]))/(w[it]*L[it]+r[it]*A[it])
        τ_r[it] = τ_w[it]
    end
    # budget-balancing social security contribution
    τ_p[it] = (pen[it]/((2+n_p[it])*(1+n_p[it-1])))/w[it]
    return nothing
end

function get_transition!(para; tol=1e-5, itermax=1000, T=25)
    K .= K[1]
    nmarket = 0
    # iterate until all markets clear
    iter = 1
    while nmarket < (T-1) && (iter <= itermax)
        #get prices, decisions, and quantities
        for it = 2:T
            factor_prices!(it, para)
        end

        for it = 2:T
            decisions!(it, para)
        end
        for it = 2:T 
            # aggregation
            quantities!(it, para)
        end
        for it = 2:T
            # government budget clearing and pension tax rates
            government!(it, para)
        end

        # check for number of markets in equilibrium
        ii = 2:T
        nmarket = sum( abs.(Y[ii]-C[ii]-I[ii]-G[ii])./Y[ii] .< tol)
        max_err = maximum(abs.(Y[ii]-C[ii]-I[ii]-G[ii])./Y[ii])
        @printf "Max error in market clearing %.04f\n" max_err
        iter += 1
    end
    if (iter >= itermax)
        @printf "No equilibrium found"
    end
end

function eqs_initial(x, para; tax=1)
    x = sqrt.(x.^2)
    K[1] = x[1]

    if tax == 1
        τ_c[1] = x[2]
    end

    # get prices, decisions, and quantities
    factor_prices!(1, para)
    decisions!(1, para)
    quantities!(1, para)

    out = zeros(2)
    out[1]  = (K[1] + B[1] - A[1])/(K[1]+1e-6)
    out[2] = (τ_c[1]*C[1] + τ_w[1]*w[1]*L[1] + τ_r[1]*r[1]*A[1] - (r[1]-n_p[1])*B[1] - G[1])/(C[1]+1e-6)
    return out
end


#function market_clearing(x, para; tax=1, T=25, tol=1e-5, damp=0.25)
para = Para()
@unpack γ, egam, β, α, δ, g = para
# Guess over capital stock and taxes
# One tax must endogenously satisfy government budget constraint
#x .= (x.^2).^(1/2)
# policies


tax .= ones(Int, T)
tax[2:end] .= 2
# solve initial equilibrium
L .= (2 .+n_p)./(1 .+n_p)
τ_p = @. κ/((2+n_p)*(1+n_p))
x = zeros(2)
x .= 0.3

f(z) = eqs_initial(z, para)
sol = nlsolve(f, x, method=:anderson)
x = sqrt.(sol.zero.^2)
println("Convergence = $(sol.f_converged)")

get_transition!(para, T=25)
@btime get_transition!(para, T=25)
# Utility

I[1] = (n_p[1]+δ)*K[1]
diff = Y - C - G - I
summ = [Y C G I diff τ_w τ_c τ_r c' K w r]
summ_dat = DataFrame(summ, :auto)
summ_dat .= round.(summ_dat, digits=2)
rename!(summ_dat, [:Y, :C, :G, :I, :diff, :τ_w, :τ_c, :τ_r, :c1, :c2, :c3, :K, :w, :r])
print(summ_dat)


para = Para()
x = zeros(2)
x .= 0.7
summ_dat = market_clearing(x, para)

