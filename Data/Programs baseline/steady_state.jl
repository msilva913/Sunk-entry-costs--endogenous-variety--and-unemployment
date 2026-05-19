# =============================================================================
# Steady-State and Calibration — Gabrovski-Silva Framework
# =============================================================================
# Canonical file. 
#
# Architecture:
#   calibrate_shares(targets) → NamedTuple of parameters (ParaCalib-compatible)
#   steady_state(para)        → NamedTuple of all SS objects
#
# Key conventions:
#   - s (worker separation conditional on firm survival) is a derived parameter.
#     It is computed inside calibrate_shares as s = (τ - δ_e) / (1 - δ_e) and
#     carried in ParaCalib for downstream use. Never set s manually without
#     ensuring consistency with τ and δ_e.
#   - f, q targets are gross rates (JOLTS-measurable). Inside calibration they
#     are divided by (1 - δ_e) to recover the effective conditional rates used
#     in model equations: (1-δ_e)*f_corr = f_target.
#   - x_m is not a primitive — it is pinned by the entry cost distribution and
#     the free-entry condition. It is returned by calibrate_shares.
# =============================================================================

using Parameters, CSV, StatsBase, Statistics, Random
using NLsolve, DataFrames
using Roots, Optim, LeastSquaresOptim
using PrettyPrinting, LaTeXStrings

# =============================================================================
# Utility Functions
# =============================================================================

"""
    compute_markup(ε)
Markup μ = ε/(ε-1) from elasticity of substitution.
"""
compute_markup(ε) = ε / (ε - 1)

"""
    compute_separation_rate(δ_e, s)
Total worker separation rate: τ = 1 - (1-δ_e)(1-s).
"""
compute_separation_rate(δ_e, s) = 1 - (1 - δ_e) * (1 - s)

"""
    compute_unemployment(τ, f, δ_e)
Steady-state unemployment from flow balance.
"""
compute_unemployment(τ, f, δ_e) = τ / (τ + (1 - δ_e) * f)

# =============================================================================
# Matching Functions
# =============================================================================

"""
    jf(θ, A, η_L)
Job-finding probability p(θ) = A*θ^(1-η_L), capped at 1.
"""
jf(θ, A, η_L) = min(A * θ^(1 - η_L), 1.0)

"""
    vf(θ, A, η_L)
Vacancy-filling probability q(θ) = A*θ^(-η_L), capped at 1.
"""
vf(θ, A, η_L) = min(A * θ^(-η_L), 1.0)

"""
    θ_invert(q, A, η_L)
Recover tightness θ from vacancy-filling rate q.
"""
function θ_invert(q, A, η_L)
    q <= 0 || q > A && throw(ArgumentError("q must be in (0, A], got q=$q, A=$A"))
    return (q / A)^(-1 / η_L)
end

# =============================================================================
# Parameter Struct
# =============================================================================

"""
    ParaCalib

Model parameters. All fields should be populated via `calibrate_shares`;
do not set `s` manually (it must satisfy s = (τ - δ_e)/(1 - δ_e) at the
calibrated steady state).
"""
@with_kw mutable struct ParaCalib
    f_e::Float64  = 1.0          # Sunk entry cost
    δ::Float64    = 0.005        # Exogenous product destruction rate
    s::Float64    = 0.0          # Worker separation conditional on firm survival
                                 # (derived: s = (τ - δ_e)/(1 - δ_e); set by calibrate_shares)
    z::Float64    = 1.0          # Technology level
    b::Float64    = 0.71         # Unemployment benefit
    ϕ::Float64    = 0.5          # Worker bargaining power
    r::Float64    = 0.04 / 12    # Monthly discount rate
    σ::Float64    = 1.0          # Inverse IES (σ=1 → log utility)
    ε::Float64    = 4.0          # Elasticity of substitution across varieties
    A::Float64    = 0.5631       # Matching function scale parameter
    η_L::Float64  = 0.6          # Matching elasticity w.r.t. unemployment
    κ::Float64    = 0.2          # Per-period matching cost (flow, paid by firm)
    ξ_inv::Float64 = 1.0         # Inverse elasticity of entry cost to vacancy value
    x_m::Float64  = 113.16       # Upper bound of sunk entry cost distribution
    ψ::Float64    = 1.5          # Shape parameter of continuation cost power law
    f_m::Float64  = 10.0         # Scale parameter of continuation cost distribution
    p_0::Float64  = 0.5          # Mass weight on continuous (power law) part of F
end

# =============================================================================
# Cost Distribution
# =============================================================================

"""
    F(x, para)

Continuation cost CDF with point mass at 0:
  F(x) = (1 - p_0) + p_0*(x/f_m)^ψ   for x ∈ [0, f_m]
  F(x) = 1                             for x > f_m
"""
function F(x, para)
    @unpack p_0, f_m, ψ = para
    (f_m <= 0 || ψ <= 0) && throw(ArgumentError("f_m and ψ must be > 0"))
    (p_0 < 0 || p_0 > 1) && throw(ArgumentError("p_0 must be in [0,1]"))
    return x < f_m ? 1 - p_0 + p_0 * (x / f_m)^ψ : 1.0
end

"""
    F_inv(F_val, f_m, ψ, p_0)
Inverse CDF accounting for atom at 0 with probability (1 - p_0).
"""
function F_inv(F_val, f_m, ψ, p_0)
    F_val <= 1 - p_0 && return 0.0
    ζ = (F_val - (1 - p_0)) / p_0
    return f_m * ζ^(1 / ψ)
end

"""
    F_inv_simple(F_val, f_m, ψ)
Inverse CDF for pure power law (no atom): x = f_m * F_val^(1/ψ).
"""
function F_inv_simple(F_val, f_m, ψ)
    (f_m <= 0 || ψ <= 0) && throw(ArgumentError("f_m and ψ must be > 0"))
    return F_val^(1 / ψ) * f_m
end

"""
    δ_e_fun(x_c, para)
Endogenous destruction rate given cost cutoff: δ_e = 1 - (1-δ)*F(x_c).
"""
δ_e_fun(x_c, para) = 1 - (1 - para.δ) * F(x_c, para)

"""
    x_c_dest_fun(δ_e, para)
Cost cutoff consistent with destruction rate δ_e (inverse of δ_e_fun).
"""
function x_c_dest_fun(δ_e, para)
    @unpack δ, f_m, ψ, p_0 = para
    F_target = (1 - δ_e) / (1 - δ)
    F_target <= 1 - p_0 && return 0.0
    ζ = (F_target - (1 - p_0)) / p_0
    return f_m * ζ^(1 / ψ)
end

"""
    dest_elast(para, x_c)
Elasticity of δ_e with respect to x_c: ε_dest = (∂δ_e/∂x_c)*(x_c/δ_e).
Returns 0 at the upper bound of the distribution.
"""
function dest_elast(para, x_c)
    @unpack ψ, f_m = para
    ζ = min((x_c / f_m)^ψ, 1.0)
    ζ >= 1.0 && return 0.0
    return ψ * ζ / (1 - ζ)
end

# =============================================================================
# Equilibrium Functions
# =============================================================================

"""
    L_fun(θ, δ_e, para)
Steady-state employment from flow balance.
"""
function L_fun(θ, δ_e, para)
    @unpack s, A, η_L = para
    f = jf(θ, A, η_L)
    τ = compute_separation_rate(δ_e, s)
    return (1 - δ_e) * f / (τ + (1 - δ_e) * f)
end

"""
    e_fun(θ, δ_e, para)
Entry rate (new vacancies posted per unit time).
"""
function e_fun(θ, δ_e, para)
    @unpack s, A, η_L = para
    f = jf(θ, A, η_L)
    τ = compute_separation_rate(δ_e, s)
    return δ_e * (θ * τ + (1 - δ_e) * f) / (τ + (1 - δ_e) * f)
end

"""
    K_fun(θ, δ_e, para)
Value of a posted vacancy (expected discounted return to posting).
"""
function K_fun(θ, δ_e, para)
    @unpack r, x_m, ξ_inv = para
    e = e_fun(θ, δ_e, para)
    Q = e^ξ_inv * x_m
    return (r + δ_e) / (1 + r) * Q
end

"""
    profit_share(δ_e, para)
Steady-state retail profit share π_s consistent with cost distribution.
"""
function profit_share(δ_e, para)
    @unpack r, ψ, p_0 = para
    μ = compute_markup(para.ε)
    cons = (ψ / (1 + ψ)) * p_0
    return (μ - 1) / μ * (1 - cons) * (r + δ_e) / (r + δ_e + cons * (1 - δ_e))
end

"""
    w_fun(θ, N, δ_e, para)
Nash wage given tightness θ, firm mass N, destruction rate δ_e.
"""
function w_fun(θ, N, δ_e, para)
    @unpack ϕ, b, z, A, η_L, κ = para
    q   = vf(θ, A, η_L)
    μ   = compute_markup(para.ε)
    w_int = N^(1 / (para.ε - 1)) * z / μ
    K   = K_fun(θ, δ_e, para)
    return (1 - ϕ) * b + ϕ * (w_int - K + θ * (K + q * κ))
end

# =============================================================================
# Equilibrium Curves (for visualization / diagnostics)
# =============================================================================

"""
    N_jcc(θ, δ_e, para)
Job creation curve: firm mass N consistent with wage bargaining at (θ, δ_e).
"""
function N_jcc(θ, δ_e, para)
    @unpack r, A, η_L, ϕ, z, b, κ, s = para
    τ   = compute_separation_rate(δ_e, s)
    q   = vf(θ, A, η_L)
    K   = K_fun(θ, δ_e, para)
    μ   = compute_markup(para.ε)
    w_int = (1 / ((1 - δ_e) * (1 - ϕ))) * (κ * q + K) * ((r + τ) / q + (1 - δ_e) * ϕ * θ) + K + b
    ρ   = w_int * (μ / z)
    return ρ^(para.ε - 1)
end

"""
    N_res(θ, δ_e, para)
Resource constraint curve: firm mass N from free entry and profit conditions.
"""
function N_res(θ, δ_e, para)
    @unpack r, f_e, ε, z, s, A, η_L = para
    τ   = compute_separation_rate(δ_e, s)
    f   = jf(θ, A, η_L)
    u   = compute_unemployment(τ, f, δ_e)
    L   = 1 - u
    π_s = profit_share(δ_e, para)
    μ   = compute_markup(ε)
    return z * L * (1 - δ_e) / (f_e * (δ_e + (r + δ_e) / (μ * π_s)))
end

"""
    x_c_exit_thresh_fun(δ_e, ρ, para)
Cost cutoff from joint exit + free-entry condition (function of ρ and δ_e).
"""
function x_c_exit_thresh_fun(δ_e, ρ, para)
    @unpack f_e, r = para
    μ   = compute_markup(para.ε)
    π_s = profit_share(δ_e, para)
    return (ρ * f_e / μ) * (1 + (r + δ_e) / (para.ε * π_s * (1 - δ_e)))
end

"""
    solve_δ_e(ρ, para)
Solve for δ_e satisfying x_c_dest = x_c_exit_thresh at given ρ.
"""
function solve_δ_e(ρ, para)
    function loss(δ_e)
        x1 = x_c_dest_fun(δ_e, para)
        x2 = x_c_exit_thresh_fun(δ_e, ρ, para)
        return (x1 - x2) / (x1 + x2)
    end
    return find_zero(loss, (1e-6, 0.3))
end

# =============================================================================
# Steady-State Solver
# =============================================================================

"""
    steady_state(para; init=0.51)

Solve for the steady state given parameters `para`.

Solves a 2×2 system in (log θ, log x_c):
  1. Job creation condition (wage bargaining)
  2. Exit condition (cost cutoff consistency)

Returns a NamedTuple of all steady-state objects.
"""
function steady_state(para; init=0.51)
    @unpack f_e, δ, s, z, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, x_m, f_m, ψ, p_0 = para

    μ    = compute_markup(ε)
    cons = (ψ / (1 + ψ)) * p_0

    function loss(y)
        θ   = exp(y[1])
        x_c = exp(y[2])

        δ_e  = δ_e_fun(x_c, para)
        τ    = compute_separation_rate(δ_e, s)
        # Profit share π_s 
        π_s  = (μ - 1) / μ * (1 - cons) * (r + δ_e) / (r + δ_e + cons * (1 - δ_e))

        f    = jf(θ, A, η_L)
        q    = vf(θ, A, η_L)
        K    = K_fun(θ, δ_e, para)
        L    = L_fun(θ, δ_e, para)
        u    = 1 - L

        # Job Creation Condition (JCC)
        lhs_jcc = (κ + K / q) * (r + τ + (1 - δ_e) * ϕ * q * θ)

        N     = π_s * z * L * (1 - δ_e) / (f_e * ((r + δ_e) / μ + δ_e * π_s))
        ρ     = N^(1 / (ε - 1))
        w_int = ρ * z / μ
        rhs_jcc = (1 - δ_e) * (1 - ϕ) * (w_int - K - b)

        # Exit Condition
        L_c     = (r + δ_e) * L / (r + δ_e + δ_e * π_s * μ)
        Y_c     = ρ * z * L_c
        x_c_new = Y_c / (ε * N) + ρ * f_e / μ

        r1 = rhs_jcc - lhs_jcc
        r2 = log(x_c_new / x_c)

        vars = (; δ_e, f, q, θ, u, ρ, w_int, K, N, L_c, Y_c, x_c, x_c_new, π_s)
        return [r1, r2], vars
    end

    sol = LeastSquaresOptim.optimize(x -> loss(x)[1], [log(init), log(δ / 2.0)], Dogleg())
    sol.converged || @warn "Solver did not fully converge (ssr=$(sol.ssr))"
    println("Steady-state solver: converged=$(sol.converged) in $(sol.iterations) iterations")

    var = loss(sol.minimizer)[2]
    @unpack δ_e, f, q, θ, u, ρ, w_int, K, N, L_c, Y_c, x_c, π_s = var
    L  = 1 - u
    L_e = L - L_c
    N_e = δ_e * N / (1 - δ_e)

    # Labor market aggregates
    v      = θ * u
    e      = δ_e * (v + 1 - u)
    Q      = e^ξ_inv * x_m
    X_v    = e / (1 + ξ_inv) * Q
    v_pret = v - e

    # Firm values
    ν_f = ρ * f_e / μ
    d_f = (r + δ_e) / (1 - δ_e) * ν_f
    X_c = N * cons * x_c

    # Wages and recruiting costs
    w   = ϕ * (w_int - K + θ * (K + q * κ)) + (1 - ϕ) * b
    X   = X_v + κ * v * q

    # Output
    C   = Y_c - X - X_c
    Y   = C + ν_f * N_e

    # Wealth and asset values
    J   = Q + (1 + r) / (1 - δ_e) * K / q
    M   = Q * v + J * L + (N + N_e) * ν_f

    # Shares
    labor_prod         = Y / (ρ * L) # remove variety effects
    labor_share        = w * L / Y
    cons_share         = C / Y
    inv_new_firm_share = ν_f * N_e / Y
    vacancy_share      = X / Y
    fixed_cost_share   = X_c / Y
    sunk_vac_cost_share = X_v / Y
    entrant_vac_share  = e / v

    x_v             = (K / q) / (κ + K / q)
    search_wedge    = w / w_int
    recruiter_share = w_int * L / Y
    ann_int_rate    = (1 + r)^12 - 1

    profit_share_ret = π_s * Y_c / Y
    profit_share_rec = ((w_int - w) * L - X) / Y

    dest_end_frac = (δ_e - δ) / δ_e
    dest_el       = dest_elast(para, x_c)

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
# Calibration
# =============================================================================

"""
Default calibration targets. Override individual fields via named-tuple merge:
    my_targets = (TARGETS..., Xc_Y=0.07)

Parameter taxonomy (must match draft Section 5.2):

  HARD-CODED TARGETS (Θ_d) — pinned by data moments, not estimated:
    X_Y:      Recruiting cost share of GDP = 1.5%. Hard-coded target in Stage 4; pins Q.
              Christiano, Eichenbaum & Trabandt (2016) estimate total hiring expenditure
              at 1–2% of GDP. 1.5% is pure vacancy/recruiting expenditure (excluding
              training). This is a calibration target, not an SMM moment.
    Xc_Y:     Fixed cost share of GDP = 10%. Hard-coded target in Stage 2; pins π_s.
              Abraham, Bormans, Konings & Roeger (2019): fixed costs ~10–15% of EU value
              added. U.S. overhead labor (Bils & Klenow 2004): ~6% of manufacturing output.
              10% is an upper bound — keep conservatively low to avoid profit squeeze.
              WARNING: π_s = 1/ε − Xc/Yc. With ε=4.3, feasible range is Xc_Y ≲ 0.13.
              Sensitivity: run (TARGETS..., Xc_Y=0.07) for lower bound.
    dest_ann: Annual product destruction rate — BED Deaths, employment-weighted.
    f, q:     Gross JOLTS rates; corrected by (1-δ_e) inside Stage 1 calibration.
    sep:      Aggregate separation rate — Shimer (2012) / JOLTS.
    η_L:      Matching elasticity — Petrongolo & Pissarides (2001), fixed externally.
    r_ann:    Annual discount rate — fixed externally.
    ε:        Elasticity of substitution — BGM / Compustat, fixed externally.
    N, w:     Normalizations.

  ESTIMATED (Θ_e) — drawn from SMM posterior; passed in as arguments, not TARGETS:
    dest_end_frac (ω_δ), p_0, b_ratio, x_v, ξ_inv, σ, and shock processes ρ_x, σ_x.
    The TARGETS tuple holds prior/default values for these during steady-state checks
    and standalone calibration runs; they are overridden by the SMM sampler at runtime.

  NOTE on p_0 and ω_δ identification:
    - p_0 (mass on continuous part of F) is not pinned by any long-run aggregate moment.
      Within Stage 3, only cons = (ψ/(1+ψ))*p_0 is identified given π_s. Estimating p_0
      from business-cycle moments via SMM allows ψ (the shape of the continuous part) to
      adjust independently, giving the distribution two free parameters.
    - ω_δ = (δ_e - δ)/δ_e is the endogenous share of total destruction. BED Deaths pins
      total δ_e but cannot separately identify how much is endogenous selection vs.
      exogenous closure — that split requires business-cycle comovement information.
"""
const TARGETS = (
    # --- Hard-coded calibration targets ---
    X_Y           = 0.015,  # Recruiting cost share of GDP [CET 2016] — Stage 4 target
    Xc_Y          = 0.10,   # Fixed cost share of GDP [Abraham et al. 2019] — Stage 2 target
    dest_ann      = 0.0754, # Annual product destruction rate [BED Deaths, emp-weighted]
    f             = 0.41,   # Gross job-finding rate [JOLTS] — corrected ÷(1-δ_e) in Stage 1
    η_L           = 0.6,    # Matching elasticity [Petrongolo & Pissarides 2001]
    q             = 0.8,    # Gross vacancy-filling rate [JOLTS] — corrected ÷(1-δ_e) in Stage 1
    sep           = 0.031,  # Aggregate separation rate [Shimer 2005 / JOLTS]
    r_ann         = 0.04,   # Annual discount rate
    ε             = 4.3,    # Elasticity of substitution [BGM / Compustat]
    N             = 1.0,    # Normalization: SS firm mass
    w             = 1.0,    # Normalization: SS wage
    # --- Estimated (Θ_e) — default/prior values; overridden by SMM sampler ---
    dest_end_frac = 0.5,    # Endogenous share of δ_e [ESTIMATED]
    p_0           = 0.5,    # Mass on continuous cost distribution [ESTIMATED]
    b_ratio       = 0.71,   # Replacement ratio b/w [ESTIMATED]
    x_v           = 1.0,    # Sunk share of total recruiting cost [ESTIMATED]
    ξ_inv         = 1,      # Inverse entry elasticity [ESTIMATED]
    σ             = 1.0,    # Inverse IES [ESTIMATED]
)

"""
    calibrate_shares(targets) → NamedTuple

Calibrate model parameters to match empirical targets in four sequential stages:

  Stage 1 — Direct conversions (targets → model primitives):
    dest_ann, dest_end_frac → δ_e, δ; r_ann → r; f, q (gross) → f_corr, q_corr
    (note: f_corr = f/(1-δ_e) because JOLTS rates are gross; model uses conditional rates)

  Stage 2 — Root-find Xc_Yc (fixed cost / consumption output):
    targets Xc_Y → pins π_s (retail profit share)

  Stage 3 — Root-find cons (distribution mass fraction):
    targets π_s → pins cons, then ψ, f_m

  Stage 4 — Root-find Q (sunk entry value):
    targets X_Y, x_v → pins Q, K, κ, z, f_e, x_c, x_m, ϕ

Returns 17 parameters suitable for constructing ParaCalib.
"""
function calibrate_shares(targets)
    @unpack X_Y, Xc_Y, dest_ann, dest_end_frac, p_0, f, η_L, q, sep, b_ratio,
            x_v, ξ_inv, ε, r_ann, σ, N, w = targets

    # ------------------------------------------------------------------
    # Stage 1: Direct conversions
    # ------------------------------------------------------------------
    δ_e      = 1 - (1 - dest_ann)^(1 / 12)
    δ        = (1 - dest_end_frac) * δ_e
    surv_prob = (1 - δ_e) / (1 - δ)   # = F(x_c) at SS

    μ  = compute_markup(ε)
    τ  = sep
    # Worker separation conditional on firm survival; derived, not a free parameter
    s  = (τ - δ_e) / (1 - δ_e)
    r  = (1 + r_ann)^(1 / 12) - 1

    # Correct gross JOLTS rates to conditional-on-survival rates used in the model:
    #   (1 - δ_e) * f_corr = f_target  ⟹  f_corr = f_target / (1 - δ_e)
    f_corr = f / (1 - δ_e)
    q_corr = q / (1 - δ_e)
    θ  = f_corr / q_corr
    u  = τ / (τ + (1 - δ_e) * f_corr)
    v  = θ * u
    L  = 1 - u
    A  = f_corr / θ^(1 - η_L)

    e   = δ_e * (v + 1 - u)
    N_e = δ_e / (1 - δ_e) * N
    b   = b_ratio * w
    ρ   = N^(1 / (ε - 1))

    # ------------------------------------------------------------------
    # Stage 2: Solve for profit share via fixed-cost ratio
    # ------------------------------------------------------------------
    # Relationship: Xc_Y = Xc_Yc * Yc_YG * (1/YG_Y)
    #   Xc_Yc ≡ X_c/Y_c (endogenous fixed cost share within consumption sector)
    #   Yc_YG = (r+δ_e)/(r+δ_e+δ_e*π_s)
    #   YG_Y  = 1 + X_Y + Xc_Y (gross output / GDP)
    function loss_xc(xc_yc)
        π_s_trial = 1 / ε - xc_yc
        Yc_YG = (r + δ_e) / (r + δ_e + δ_e * π_s_trial)
        YG_Y  = 1 + X_Y + Xc_Y
        return Xc_Y / (Yc_YG * YG_Y) - xc_yc
    end
    Xc_Yc = find_zero(loss_xc, (0.01, 0.5))
    π_s   = 1 / ε - Xc_Yc

    L_c = (r + δ_e) * L / (r + δ_e + δ_e * π_s * μ)
    L_e = L - L_c

    # ------------------------------------------------------------------
    # Stage 3: Solve for cost distribution mass fraction
    # ------------------------------------------------------------------
    # π_s(cons) = [(μ-1)/μ] * (1-cons) * (r+δ_e) / (r+δ_e + cons*(1-δ_e))
    function loss_psi(cons_trial)
        return (μ - 1) / μ * (1 - cons_trial) * (r + δ_e) /
               (r + δ_e + cons_trial * (1 - δ_e)) - π_s
    end
    cons = find_zero(loss_psi, 0.1)
    ψ_c  = cons / p_0
    ψ    = ψ_c / (1 - ψ_c)

    # Pre-compute surplus ratio for Stage 4
    surplus_ratio = (r + τ) / (1 - δ_e) / (q_corr * x_v)

    # ------------------------------------------------------------------
    # Stage 4: Solve for vacancy value Q
    # ------------------------------------------------------------------
    # Output consistency: Y implied by X/X_Y should equal Y from national accounts
    function loss_Q(Q_trial)
        Q_trial = abs(Q_trial)
        K     = Q_trial * (r + δ_e) / (1 + r)
        κ     = (1 - x_v) / x_v * K / q_corr
        X     = e / (1 + ξ_inv) * Q_trial + κ * q_corr * v   # sunk + flow recruiting
        w_int = surplus_ratio * K + w + K
        z     = (μ / ρ) * w_int
        f_e   = π_s * z * L_c * (1 - δ_e) * μ / (N * (r + δ_e))
        ν_f   = ρ * f_e / μ
        d_f   = (r + δ_e) / (1 - δ_e) * ν_f
        Y_c   = ρ * z * L_c
        x_c   = Y_c / (ε * N) + ν_f
        X_c   = N * cons * x_c
        C     = Y_c - X - X_c
        Y_new = C + ν_f * N_e
        Y_rec = X / X_Y                   # Y implied by recruiting share target
        loss_val = 100 * (Y_rec - Y_new) / (Y_rec + Y_new)
        return loss_val, (; w_int, κ, z, f_e, K, d_f, ν_f, x_c, X_c, X, C, Y_c, Y=Y_new)
    end

    Q_opt = find_zero(Q -> loss_Q(Q)[1], 1.0)
    Q     = abs(Q_opt)
    _, Q_out = loss_Q(Q)
    @unpack w_int, κ, z, f_e, K, d_f, ν_f, x_c, X_c, X, C, Y_c = Q_out

    # ------------------------------------------------------------------
    # Stage 5: Recover remaining parameters
    # ------------------------------------------------------------------
    ζ       = (surv_prob - (1 - p_0)) / p_0   # ζ = (x_c/f_m)^ψ at SS
    dest_el = ψ * ζ / (1 - ζ)
    f_m     = x_c / ζ^(1 / ψ)

    ϕ   = (w - b) / (w_int - K + θ * (K + q_corr * κ) - b)
    x_m = Q / e^ξ_inv

    return (;
        f_e, δ, z, b, ϕ, r, σ, ε, A, η_L,
        κ, ξ_inv, x_m, s, ψ, p_0, f_m
    )
end
