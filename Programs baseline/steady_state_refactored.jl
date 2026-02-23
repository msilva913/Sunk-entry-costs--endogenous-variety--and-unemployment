# =============================================================================
# Steady-State and Calibration Analysis - REFACTORED
# =============================================================================
# This refactored version improves the architecture by:
# - Consolidating repeated calculations
# - Adding comprehensive error handling
# - Separating module code from executable code
# - Extracting common patterns into utility functions
# =============================================================================

using Parameters, CSV, StatsBase, Statistics, Random
using NLsolve, DataFrames
using Roots, Optim, LeastSquaresOptim
using PrettyPrinting, LaTeXStrings

# =============================================================================
# Utility Functions - Common Calculations
# =============================================================================

"""
    compute_markup(ε)
Compute markup (μ) from elasticity of substitution.
"""
compute_markup(ε) = ε / (ε - 1)

"""
    compute_separation_rate(δ_e, s)
Compute total separation rate τ combining exogenous and endogenous destruction.
"""
compute_separation_rate(δ_e, s) = 1 - (1 - δ_e) * (1 - s)

"""
    compute_unemployment(τ, f, δ_e)
Compute unemployment rate from separation rate and job-finding rate.
"""
compute_unemployment(τ, f, δ_e) = τ / (τ + (1 - δ_e) * f)


# =============================================================================
# Matching Functions
# =============================================================================

"""
    jf(θ, A, η_L)

Job finding probability as a function of market tightness θ.
Capped at 1.0 for feasibility.
"""
jf(θ, A, η_L) = min(A * θ^(1 - η_L), 1.0)

"""
    vf(θ, A, η_L)

Vacancy filling probability as a function of market tightness θ.
Capped at 1.0 for feasibility.
"""
vf(θ, A, η_L) = min(A * θ^(-η_L), 1.0)

"""
    θ_invert(q, A, η_L)

Invert vacancy filling probability to obtain market tightness θ.
Validates inputs to ensure q ∈ (0, A].
"""
function θ_invert(q, A, η_L)
    if q <= 0 || q > A
        throw(ArgumentError("Vacancy filling rate q must be in (0, A], got q=$q, A=$A"))
    end
    return (q / A)^(-1 / η_L)
end

# =============================================================================
# Parameter Struct
# =============================================================================

@with_kw mutable struct ParaCalib
    f_e::Float64 = 1.0                 # Entry cost
    τ::Float64 = 0.0305                # Worker separation rate (exogenous)
    δ::Float64 = 0.005                 # Product destruction rate (exogenous)
    s::Float64 = 0.0                   # Separation rate adjustment
    z::Float64 = 1.0                   # Technology level
    b::Float64 = 0.71                  # Unemployment insurance
    ϕ::Float64 = 0.5                   # Bargaining power 
    r::Float64 = 0.04 / 12             # Rate of time preference (monthly)
    σ::Float64 = 1.0                   # Inverse intertemporal elasticity of substitution
    ε::Float64 = 4.0                   # Elasticity of substitution
    A::Float64 = 0.5631                # Job matching level parameter
    η_L::Float64 = 0.6                 # Elasticity of matching function w.r.t. unemployment
    κ::Float64 = 0.2                   # Fixed matching cost
    ξ_inv::Float64 = 1.0               # Inverse elasticity of entry to vacancy value
    x_m::Float64 = 113.16              # Upper bound on cost distribution
    
    # Parameters related to firm heterogeneity and cost distribution
    ψ::Float64 = 1.5                   # Curvature (shape parameter of power law)
    f_m::Float64 = 10.0                # Max value of cost draw (location parameter)
    p_0::Float64 = 0.5                 # Probability of continuous (power law) part
end

# =============================================================================
# Equilibrium Functions - Cost Distribution
# =============================================================================

"""
    F(x, para)

Continuation cost-draw cdf with atom at 0.
F(x) = 1 - p_0 + p_0*(x/f_m)^ψ for x ∈ [0, f_m]
F(x) = 1.0 for x > f_m
"""
function F(x, para)
    @unpack p_0, f_m, ψ = para
    
    # Input validation
    if f_m <= 0 || ψ <= 0
        throw(ArgumentError("f_m and ψ must be > 0, got f_m=$f_m, ψ=$ψ"))
    end
    if p_0 < 0 || p_0 > 1
        throw(ArgumentError("p_0 must be in [0,1], got p_0=$p_0"))
    end
    
    if x < f_m
        return 1 - p_0 + p_0 * (x / f_m)^ψ
    else
        return 1.0
    end
end 

"""
    F_inv_simple(F_val, f_m, ψ)

Inverse cdf for simple power law (without atom).
Solves F_val = (x/f_m)^ψ for x.
"""
function F_inv_simple(F_val, f_m, ψ)
    if f_m <= 0 || ψ <= 0
        throw(ArgumentError("f_m and ψ must be > 0, got f_m=$f_m, ψ=$ψ"))
    end
    return F_val^(1/ψ) * f_m
end

"""
    F_inv(F_val, f_m, ψ, p_0)

Inverse cdf for bounded power law with atom at 0.
Accounts for atom at 0 with probability (1-p_0).
"""
function F_inv(F_val, f_m, ψ, p_0)
    if F_val <= 1 - p_0
        return 0.0
    else
        ζ = (F_val - (1 - p_0)) / p_0
        return f_m * ζ^(1/ψ)
    end
end

"""
    δ_e_fun(x_c, para)

Endogenous product destruction rate given cutoff cost x_c.
δ_e = 1 - (1-δ)*F(x_c)
"""
function δ_e_fun(x_c, para)
    @unpack δ = para
    return 1 - (1 - δ) * F(x_c, para)
end

"""
    x_c_dest_fun(δ_e, para)

Exit threshold that achieves target destruction rate δ_e.
Solves for x_c such that δ_e = 1 - (1-δ)*F(x_c).
"""
function x_c_dest_fun(δ_e, para)
    @unpack δ, f_m, ψ, p_0 = para
    F_target = (1 - δ_e) / (1 - δ)
    if F_target <= 1 - p_0
        return 0.0
    else
        ζ = (F_target - (1 - p_0)) / p_0
        return f_m * ζ^(1/ψ)
    end
end

"""
    dest_elast(para, x_c)

Destruction rate elasticity with respect to exit threshold.
ε_dest = (∂δ_e/∂x_c) * (x_c/δ_e)
"""
function dest_elast(para, x_c)
    @unpack ψ, f_m = para 
    ζ = (x_c / f_m)^ψ
    ζ = min(ζ, 1.0)
    if ζ >= 1.0
        return 0.0  # Elasticity is zero when at upper bound
    else
        return ψ * ζ / (1 - ζ)
    end
end

# =============================================================================
# Equilibrium Functions - Labor Market
# =============================================================================

"""
    L_fun(θ, δ_e, para)

Steady-state employment as a function of market tightness θ.
"""
function L_fun(θ, δ_e, para)
    @unpack s, A, η_L = para
    f = jf(θ, A, η_L)
    τ = compute_separation_rate(δ_e, s)
    return (1 - δ_e) * f / (τ + (1 - δ_e) * f)
end

"""
    e_fun(θ, δ_e, para)

New vacancy/entry rate as a function of market tightness θ.
"""
function e_fun(θ, δ_e, para)
    @unpack s, A, η_L = para
    f = jf(θ, A, η_L)
    τ = compute_separation_rate(δ_e, s)
    return δ_e * (θ * τ + (1 - δ_e) * f) / (τ + (1 - δ_e) * f)
end

"""
    K_fun(θ, δ_e, para)

Value of a vacancy as a function of market tightness θ and destruction rate.
"""
function K_fun(θ, δ_e, para)
    @unpack r, x_m, ξ_inv = para
    e = e_fun(θ, δ_e, para)
    Q = e^ξ_inv * x_m
    return (r + δ_e) / (1 + r) * Q
end

"""
    w_fun(θ, N, δ_e, para)

Wage as a function of market tightness θ, firm mass N, and destruction rate.
"""
function w_fun(θ, N, δ_e, para)
    @unpack ϕ, b, z, A, η_L, κ = para
    q = vf(θ, A, η_L)
    μ = compute_markup(para.ε)
    w_int = N^(1 / (para.ε - 1)) * z / μ
    K = K_fun(θ, δ_e, para)
    return (1 - ϕ) * b + ϕ * (w_int - K + θ * (K + q * κ))
end

"""
    profit_share(δ_e, para)

Retail profit share consistent with cost structure.
"""
function profit_share(δ_e, para)
    @unpack r, ψ, p_0 = para
    μ = compute_markup(para.ε)
    ψ_c = ψ / (1 + ψ)
    cons = ψ_c * p_0
    π_s = (μ - 1) / μ * (1 - cons) * (r + δ_e) / (r + δ_e + cons * (1 - δ_e))
    return π_s 
end

# =============================================================================
# Equilibrium Functions - Job Creation Curves
# =============================================================================

"""
    N_jcc(θ, δ_e, para)

Job creation curve: Firm mass N as a function of θ.
Derived from wage-setting condition.
"""
function N_jcc(θ, δ_e, para)
    @unpack r, A, η_L, δ, ϕ, z, b, κ, s = para
    τ = compute_separation_rate(δ_e, s)
    q = vf(θ, A, η_L)
    K = K_fun(θ, δ_e, para)
    μ = compute_markup(para.ε)
    w_int = (1 / ((1 - δ_e) * (1 - ϕ))) * (κ * q + K) * ((r + τ) / q + (1 - δ_e) * ϕ * θ) + K + b
    ρ = w_int * (μ / z)
    N = ρ^(para.ε - 1)
    return N
end

"""
    N_res(θ, δ_e, para)

Resource constraint curve: Firm mass N as a function of θ.
Derived from free entry and profit conditions.
"""
function N_res(θ, δ_e, para)
    @unpack r, f_e, ε, z, s, A, η_L = para
    τ = compute_separation_rate(δ_e, s)
    f = jf(θ, A, η_L)
    u = compute_unemployment(τ, f, δ_e)
    L = 1 - u
    π_s = profit_share(δ_e, para)
    μ = compute_markup(ε)
    N = z * L * (1 - δ_e) / (f_e * (δ_e + (r + δ_e) / (μ * π_s)))
    return N
end

"""
    x_c_exit_thresh_fun(δ_e, ρ, para)

Exit threshold from free entry and profit maximization.
Combines exit condition with free entry condition.
"""
function x_c_exit_thresh_fun(δ_e, ρ, para)
    @unpack f_e, r = para 
    μ = compute_markup(para.ε)
    π_s = profit_share(δ_e, para)
    return (ρ * f_e / μ) * (1 + (r + δ_e) / (para.ε * π_s * (1 - δ_e)))
end

# =============================================================================
# Steady-State Solver
# =============================================================================

"""
    steady_state(para; init=0.51)

Compute steady-state statistics given parameters.
Solves system of two equilibrium conditions:
1. Job creation condition (wage-setting)
2. Exit condition (product destruction)

Returns a NamedTuple of all key statistics.
"""
function steady_state(para; init=0.51)
    @unpack f_e, δ, s, z, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, x_m, f_m, ψ, p_0 = para
    
    # Pre-compute constants
    μ = compute_markup(ε)
    ψ_c = ψ / (1 + ψ)
    cons = ψ_c * p_0

    function loss(y) # y = [log θ, log x_c]
        θ = exp(y[1])
        x_c = exp(y[2])

        δ_e = δ_e_fun(x_c, para)
        τ = compute_separation_rate(δ_e, s)

        π_s = (μ - 1) / μ * (1 - cons) * (r + δ_e) / (r + δ_e + cons * (1 - δ_e))

        f = jf(θ, A, η_L)
        q = vf(θ, A, η_L)
        K = K_fun(θ, δ_e, para)
        L = L_fun(θ, δ_e, para)
        u = 1 - L

        # Job Creation Condition
        lhs_jcc = (κ + K / q) * (r + τ + (1 - δ_e) * ϕ * q * θ)

        N = π_s * z * L * (1 - δ_e) / (f_e * ((r + δ_e) / μ + δ_e * π_s))
        ρ = N^(1 / (ε - 1))
        w_int = ρ * z / μ
        rhs_jcc = (1 - δ_e) * (1 - ϕ) * (w_int - K - b)

        # Exit Condition
        L_c = (r + δ_e) * L / (r + δ_e + δ_e * π_s * μ)
        Y_c = ρ * z * L_c
        x_c_new = Y_c / (ε * N) + ρ * f_e / μ

        r1 = rhs_jcc - lhs_jcc
        r2 = log(x_c_new / x_c)

        resid = [r1, r2]
        vars = (; δ_e, f, q, θ, u, ρ, w_int, K, N, L_c, Y_c, x_c, x_c_new, π_s)
        return resid, vars
    end

    # Solve using least squares
    sol = LeastSquaresOptim.optimize(x -> loss(x)[1], [log(init), log(δ / 2.0)], Dogleg())
    
    if !sol.converged
        @warn "Solver did not fully converge: converged=$(sol.converged)"
    end
    
    println("Steady-state solver: converged=$(sol.converged) in $(sol.iterations) iterations")
    
    var = loss(sol.minimizer)[2]
    @unpack δ_e, f, q, θ, u, ρ, w_int, K, N, L_c, Y_c, x_c, π_s = var
    L = 1 - u
    L_e = L - L_c
    N_e = δ_e * N / (1 - δ_e)

    # Labor market variables
    v = θ * u
    e = δ_e * (v + 1 - u)
    Q = e^ξ_inv * x_m
    X_v = e * (1 / (1 + ξ_inv)) * Q
    v_pret = v - e

    # Firm values and costs
    ν_f = ρ * f_e / μ
    X_c = N * cons * x_c
    d_f = (r + δ_e) / (1 - δ_e) * ν_f

    # Wages and costs
    w = ϕ * (w_int - K + θ * (K + q * κ)) + (1 - ϕ) * b
    X = X_v + κ * v * q

    # Output and aggregate variables
    C = Y_c - X - X_c
    Y = C + ν_f * N_e

    J = Q + (1 + r) / (1 - δ_e) * K / q
    M = Q * v + J * L + (N + N_e) * ν_f
    labor_prod = Y / (ρ * L)

    # Shares
    labor_share = w * L / Y
    cons_share = C / Y
    inv_new_firm_share = ν_f * N_e / Y
    vacancy_share = X / Y
    fixed_cost_share = X_c / Y
    sunk_vac_cost_share = X_v / Y
    entrant_vac_share = e / v

    x_v = (K / q) / (κ + K / q)
    search_wedge = w / w_int
    recruiter_share = w_int * L / Y
    ann_int_rate = (1 + r)^12 - 1

    # Profit decomposition
    profit_share_ret = π_s * Y_c / Y
    profit_share_rec = ((w_int - w) * L - X) / Y

    dest_end_frac = (δ_e - δ) / δ_e
    dest_el = dest_elast(para, x_c)

    return (;
        θ, δ_e, x_c, N, f, q, u, v, v_pret, e, K, ρ, N_e,
        ν_f, d_f, w_int, w, L, L_e, L_c, Q, J,
        X_v, X, X_c, C, Y_c, Y, labor_share, dest_end_frac,
        labor_prod, cons_share, inv_new_firm_share, vacancy_share, sunk_vac_cost_share, 
        M, entrant_vac_share, x_v, search_wedge, recruiter_share, μ, ann_int_rate,
        π_s, profit_share_rec, profit_share_ret, dest_el
    )
end

# =============================================================================
# Calibration Function
# =============================================================================

# Default calibration targets
const TARGETS = (
    X_Y=0.015,            # recruiting cost share of output
    Xc_Y=0.1,             # fixed cost share of output
    dest_ann=0.0754,      # annual product destruction rate
    dest_end_frac=0.5,    # endogenous share of destruction rate
    p_0=0.5,              # probability of continuous part
    f=0.41,               # job-finding rate
    η_L=0.6,              # elasticity of matching fun w.r.t. unemployment
    q=0.8,                # vacancy filling rate
    sep=0.031,            # aggregate separation rate
    b_ratio=0.71,         # ratio of unemployment benefits to wage
    x_v=1.0, 
    ξ_inv=1,              # inverse elasticity of entry to vacancy value
    r_ann=0.04,           # annual discount rate
    ε=4.3,                # elasticity of substitution
    σ=1.0,                # inverse IES
    N=1.0,                # SS mass of firms (normalization)
    w=1.0,                # SS wage (normalization)
)

"""
    calibrate_shares(targets)

Calibrate model parameters to match empirical targets.

Calibration proceeds in three sequential stages:
  1. **Stage 1 (Direct conversions)**: Convert targets to intermediate labor-market variables
     - destruction rates: dest_ann, dest_end_frac → δ_e, δ
     - time preference: r_ann → r (monthly)
     - labor market: f, q, η_L, sep → θ, u, v, L, A
  
  2. **Stage 2 (Profit share)**: Root-find consumption-sector output share
     - Target: Xc_Y (fixed cost share)
     - Pin: Xc_Yc (fixed cost share within consumption sector), then π_s (profit share)
  
  3. **Stage 3 (Distribution parameters)**: Root-find consumption parameter
     - Target: π_s (from Stage 2)
     - Pin: cons (cost distribution mass fraction), then ψ, f_m
  
  4. **Stage 4 (Vacancy value)**: Root-find vacancy value Q
     - Targets: X_Y (recruiting share), x_v (vacancy-value split)
     - Pin: Q, then K, κ, w_int, z, f_e, x_c, x_m, ϕ

Returns a NamedTuple of 18 calibrated parameters for ParaCalib.
"""
function calibrate_shares(targets)
    @unpack X_Y, Xc_Y, dest_ann, dest_end_frac, p_0, f, η_L, q, sep, b_ratio, x_v, ξ_inv, ε, r_ann, σ, N, w = targets

    # =================================================================
    # STAGE 1: Direct Conversions (Targets → Fundamentals)
    # =================================================================
    # Target: dest_ann (annual rate) → δ_e (monthly destruction rate)
    # Target: dest_end_frac → fraction of δ_e that is endogenous
    δ_e = 1 - (1 - dest_ann)^(1 / 12)
    δ = (1.0 - dest_end_frac) * δ_e
    surv_prob = (1 - δ_e) / (1 - δ)
    
    # Target: ε → markup; sep → worker separation rate
    μ = compute_markup(ε)
    τ = sep
    s = (τ - δ_e) / (1 - δ_e)
    
    # Target: r_ann → monthly discount rate
    r = (1 + r_ann)^(1 / 12) - 1

    # Targets: f, q (corrected for destruction), η_L → labor market tightness, unemployment
    # Note: f and q are corrected rates (conditional on survival); divide out survival prob
    f_corr = f / (1 - δ_e)
    q_corr = q / (1 - δ_e)
    θ = f_corr / q_corr
    u = τ / (τ + (1 - δ_e) * f_corr)
    v = θ * u
    L = 1 - u

    # Target: η_L → matching function level
    A = f_corr / θ^(1 - η_L)

    # Derived from job creation/destruction balance
    e = δ_e * (v + 1 - u)
    N_e = δ_e / (1 - δ_e) * N
    b = b_ratio * w
    ρ = N^(1 / (ε - 1))

    # =================================================================
    # STAGE 2: Solve for Profit Share (Root-find Xc_Yc)
    # =================================================================
    # Target: Xc_Y (fixed cost as share of total output)
    # Relationship: Xc_Y = (Xc_Yc) * (Yc_YG) * (YG_Y)
    #   where Xc_Yc ≡ fixed cost / consumption output (endogenous)
    #         Yc_YG ≡ consumption output / goods output (from national accounts)
    #         YG_Y  ≡ goods output / total output (from national accounts)
    # Solving for Xc_Yc pins the profit share π_s.
    
    function loss_xc(x)
        # x is Xc_Yc; solve π_s = 1/ε - x (profit share = retail markup - fixed cost ratio)
        π_s_trial = 1 / ε - x
        # Fraction of retail output spent on consumption vs goods
        Yc_YG = (r + δ_e) / (r + δ_e + δ_e * π_s_trial)
        # Gross output is consumption + fixed costs + recruiting
        YG_Y = 1 + X_Y + Xc_Y
        # Match: target ratio should equal computed ratio
        return Xc_Y / (Yc_YG * YG_Y) - x
    end 
    
    # Bracket for search: Xc_Yc ∈ [0.01, 0.5] (fixed cost ratio must be positive, less than 50% of profit share)
    Xc_Yc = find_zero(loss_xc, (0.01, 0.5))
    π_s = 1 / ε - Xc_Yc
    
    # Allocate labor between consumption-production and goods-production
    L_c = (r + δ_e) * L / (r + δ_e + δ_e * π_s * μ)
    L_e = L - L_c

    # =================================================================
    # STAGE 3: Solve for Cost Distribution Parameter (Root-find cons)
    # =================================================================
    # Target: π_s (from Stage 2)
    # Relationship: π_s(cons) = [(μ-1)/μ] * (1-cons) * (r+δ_e) / (r+δ_e + cons*(1-δ_e))
    # Solving for cons (mass fraction of distribution with costs) pins ψ.
    
    function loss_psi(cons_trial)
        π_s_new = (μ - 1) / μ * (1 - cons_trial) * (r + δ_e) / (r + δ_e + cons_trial * (1 - δ_e))
        return π_s_new - π_s 
    end 
    
    cons = find_zero(loss_psi, 0.1)
   
    
    # Convert consumption parameter to distribution shape: ψ = cons / (1 - cons) * (1 / p_0)
    ψ_c = cons / p_0
    ψ = ψ_c / (1 - ψ_c)

    # Pre-compute wage-setting surplus ratio for later use
    surplus_ratio = (r + τ) / (1 - δ_e) * (1 / (q_corr * x_v))

    # =================================================================
    # STAGE 4: Solve for Vacancy Value (Root-find Q)
    # =================================================================
    # Targets: X_Y (recruiting cost share), x_v (share of vacancy value allocated to sunk costs)
    # Relationship: X_Y = [recruiting flows] / [total output]
    # Complex equilibrium: Q → K → w_int → z → f_e → Y via aggregate consistency.
    
    function loss_Q(Q_trial)
        Q_trial = abs(Q_trial)
        # Step A: Vacancy value → contact rate value
        K = Q_trial * (r + δ_e) / (1 + r)
        
        # Step B: Contact value split into fixed vs variable matching cost
        κ = (1 - x_v) / x_v * K / q_corr
        
        # Step C: Recruiting expenditure = sunk entry + variable matching cost
        X = e / (1 + ξ_inv) * Q_trial + κ * q_corr * v
        
        # Step D: Wage-setting interior solution
        w_int = surplus_ratio * K + w + K
        
        # Step E: Firm productivity required to support w_int
        z = (μ / ρ) * w_int
        
        # Step F: Entry cost from free-entry condition
        f_e = π_s * z * L_c * (1 - δ_e) * μ / (N * (r + δ_e))
        
        # Step G: Firm value and fixed costs
        ν_f = ρ * f_e / μ
        d_f = (r + δ_e) / (1 - δ_e) * ν_f
        Y_c = ρ * z * L_c
        x_c = Y_c / (ε * N) + ν_f
        X_c = N * cons * x_c
        
        # Step H: Output accounting: Y = C + ν_f*(N_e)
        C = Y_c - X - X_c
        Y_new = C + ν_f * N_e
        
        # Step I: Consistency check: recruiting cost share of output
        Y_from_recruiting = 1 / X_Y * X
        
        # Loss: output implied by recruiting costs should match computed output
        loss_val = 100 * (Y_from_recruiting - Y_new) / (Y_from_recruiting + Y_new)
        
        # Return loss and all intermediate values for extraction
        return loss_val, (; w_int, κ, z, f_e, K, d_f, ν_f, x_c, X_c, X, C, Y_c, Y=Y_new)
    end

    # Root-find: Q ∈ [0.001, 1.0] with robust bracketing
    loss_Q_residual(Q) = loss_Q(Q)[1]
    Q_opt = find_zero(loss_Q_residual, 1.0)
    Q = abs(Q_opt)
    
    # Extract all computed values from loss function (avoid double-evaluation)
    _, Q_results = loss_Q(Q)
    @unpack w_int, κ, z, f_e, K, d_f, ν_f, x_c, X_c, X, C, Y_c = Q_results

    # =================================================================
    # STAGE 5: Final Parameter Calculations
    # =================================================================
    # Compute distribution quantile and upper bound
    ζ = (surv_prob - (1 - p_0)) / p_0
    dest_el = ψ * ζ / (1 - ζ)
    f_m = x_c / ζ^(1 / ψ)

    # Extract bargaining power from wage equation
    w_surplus = w_int - K + θ * (K + q_corr * κ) - b  # Surplus available for wage negotiation
    ϕ = (w - b) / w_surplus
    
    # Upper bound of cost distribution
    x_m = Q / e^ξ_inv

    return (;
        f_e, δ=δ, z=z, b=b, ϕ=ϕ, r, σ=σ, ε=ε, A=A, η_L=η_L,
        κ=κ, ξ_inv=ξ_inv, x_m=x_m, s=s, ψ=ψ, p_0=p_0, f_m=f_m
    )
end

# =============================================================================
# Helper: Solver for endogenous destruction rate
# =============================================================================

"""
    solve_δ_e(ρ, para)

Solve for endogenous destruction rate δ_e that satisfies equilibrium conditions.
Matches exit threshold condition with product destruction dynamics.
"""
function solve_δ_e(ρ, para)
    function loss(δ_e)
        x_c1 = x_c_dest_fun(δ_e, para)
        x_c2 = x_c_exit_thresh_fun(δ_e, ρ, para)
        return (x_c1 - x_c2) / (x_c1 + x_c2)
    end 
    δ_e = find_zero(loss, (1e-6, 0.3))
    return δ_e
end

