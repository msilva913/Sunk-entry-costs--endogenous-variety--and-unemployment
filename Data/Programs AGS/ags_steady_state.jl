# ags_steady_state.jl
# =============================================================================
# Calibration and steady state of the AGS model with costly reactivation (Part 1 of
# decision D15; equations in Notes/ags_equilibrium.md). Standalone: no variety, no product-line
# entry block, no endogenous exit. Destroyed product lines are replaced by clones at no cost.
#
# Economy (monthly). Linear utility (σ = 0), so the discount factor is β = 1/(1+r).
#   Output per worker z̄·z. No markup: there is no variety, so z is measured productivity.
#   Matching q(θ) = Aθ^(−η), f(θ) = θ q(θ).
#   Positions are durable. Creation: Q = x_m e^ξ_inv (marginal cost of the e-th new position).
#   Per-hire cost κ ≥ 0. Destruction δ removes jobs and vacancies; separation s ends a match only.
#   Reactivation: a vacated position draws a cost c ~ F_r (mass 1 − p_r at zero, power law with
#   shape ψ_r on [0, a]) and is reactivated iff c ≤ Q. Λ_r(Q) = min{(1−p_r) + p_r (Q/a)^ψ_r, 1}.
#
# Calibration targets (the baseline's, steady_state.jl TARGETS, where they apply):
#   dest_ann → δ (= measured exit: no endogenous exit), sep → τ, f and q (gross JOLTS, divided by
#   1 − δ as in the baseline), η_L, r_ann, the wage normalization w = 1, b/w, the sunk share x_v
#   of recruiting cost, the recruiting cost share X_Y, ξ_inv, and Λ̄_r with p_r, ψ_r.
# Pinned: δ, s, A, θ, u, v, e, then Q in closed form from X_Y (z̄ and X are both linear in Q),
#   and K, κ, z̄, ϕ, x_m, a.
#
# The steady state is solved again from the parameters alone (steady_state_ags) so that
# ags_steady_state_checks.jl can test the round trip.
# =============================================================================

using Parameters, Roots, Printf

"""
Default AGS targets. Values shared with the baseline match steady_state.jl TARGETS.
Λr_target = nothing gives costless reactivation (Λ_r ≡ 1: Gabrovski–Silva with three shocks).
The reactivation defaults are the D15 working point (Λ̄_r = 0.25, p_r = 0.8, ψ_r = 1).
"""
const AGS_TARGETS = (
    dest_ann  = 0.0320,   # BED Deaths, employment-weighted, 1993–2019 (D1)
    sep       = 0.031,    # aggregate separation rate τ, monthly
    f         = 0.41,     # gross job-finding rate (JOLTS), monthly
    q         = 0.8,      # gross vacancy-filling rate (JOLTS), monthly
    η_L       = 0.6,      # matching elasticity
    r_ann     = 0.04,     # annual discount rate
    w         = 1.0,      # normalization: steady-state wage
    b_ratio   = 0.9,      # b/w
    x_v       = 0.5,      # sunk share of recruiting cost: K/q over (κ + K/q)
    X_Y       = 0.015,    # recruiting cost share of output (Part 1: the data target)
    ξ_inv     = 1.0,      # inverse elasticity of vacancy creation (0 = vacancy free entry)
    Λr_target = 0.25,     # steady-state reactivation rate; nothing = costless
    p_r       = 0.8,      # mass of the costly part of F_r
    ψ_r       = 1.0,      # shape of F_r's power law
)

"""
Recruiting-cost target for the Part 1 → Part 2 bridge: the baseline's X/(w_int·L), so AGS and
the full model share recruiting cost relative to match revenue. In AGS, output is w_int·L, so
this is AGS's X/Y. Value: 1.5 % × Y/(w_int·L) = 1.5 % × 1.1165 = 1.67 %, from the baseline
steady state (steady_state.jl, TARGETS with p_0 = 0.8; the ratio is the same with costless or
costly reposting and at ξ_inv = 1 or 0.01). Scratch computation, Nov 4, 2026: the bridge
program must recompute it from the baseline at the same draw before any number is reported.
Use: merge(AGS_TARGETS, (X_Y = X_Y_BRIDGE,)).
"""
const X_Y_BRIDGE = 0.0167

# ── Reactivation ────────────────────────────────────────────────────────────

"""Λ_r = F_r(Q) = min{(1−p_r) + p_r (Q/a)^ψ_r, 1}; 1 under costless reactivation (a = 0)."""
Λr_ags(Q, cal) = cal.a == 0 ? 1.0 : min(1 - cal.p_r + cal.p_r * (Q / cal.a)^cal.ψ_r, 1.0)

"""E max(Q − c, 0) per vacated position, interior closed form: Q[Λ_r − ψ_r/(1+ψ_r)(Λ_r − 1 + p_r)]."""
Qrep_ags(Q, Λ_r, cal) = cal.a == 0 ? Q : Q * (Λ_r - cal.ψ_r / (cal.ψ_r + 1) * (Λ_r - 1 + cal.p_r))

# ── Calibration ─────────────────────────────────────────────────────────────

"""
    calibrate_ags(tg) → NamedTuple

Parameters of the AGS model from targets `tg` (see AGS_TARGETS). Closed form; no root-finding.
"""
function calibrate_ags(tg = AGS_TARGETS)
    @unpack dest_ann, sep, f, q, η_L, r_ann, w, b_ratio, x_v, X_Y, ξ_inv, Λr_target, p_r, ψ_r = tg
    repost = !isnothing(Λr_target)
    if repost
        p_r > 0 || error("Λr_target needs p_r > 0.")
        1 - p_r < Λr_target < 1 || error("Λr_target = $Λr_target must lie inside (1−p_r, 1) = ($(1-p_r), 1).")
    end
    Λ_r = repost ? Λr_target : 1.0

    # Rates, flows and stocks
    δ  = 1 - (1 - dest_ann)^(1 / 12)
    r  = (1 + r_ann)^(1 / 12) - 1
    τ  = sep
    s  = (τ - δ) / (1 - δ)
    fc = f / (1 - δ)                       # rates conditional on survival, as in the baseline
    qc = q / (1 - δ)
    θ  = fc / qc
    u  = τ / (τ + (1 - δ) * fc)
    v  = θ * u
    L  = 1 - u
    A  = fc / θ^(1 - η_L)
    e  = δ * (v + L) + (1 - δ) * (1 - Λ_r) * s * L   # replaces destroyed and retired positions

    # Everything per unit Q: K = kQ, κ = κQ·Q, 𝓡 = gQ·Q
    k    = (r + δ) / (1 + r)
    κQ   = (1 - x_v) / x_v * k / qc
    gap  = repost ? (1 - Λ_r) + ψ_r / (ψ_r + 1) * (Λ_r - 1 + p_r) : 0.0
    gQ   = (1 - δ) * gap / (1 + r)
    # Job creation in steady state: z̄ = w + K + s𝓡 + (r+τ)/(1−δ)·(κ + K/q), linear in Q
    cw   = k + s * gQ + (r + τ) / (1 - δ) * (κQ + k / qc)
    # Recruiting costs X = eQ/(1+ξ_inv) + κ q v, linear in Q
    aX   = e / (1 + ξ_inv) + κQ * qc * v
    den  = aX - X_Y * cw * L
    den > 0 || error("X_Y = $X_Y is infeasible: recruiting costs cannot reach that share of output.")
    Q    = X_Y * w * L / den

    K    = k * Q
    κ    = κQ * Q
    R_sh = gQ * Q
    zbar = w + cw * Q
    b    = b_ratio * w
    ϕ    = (w - b) / (zbar - K + θ * (K + qc * κ) - s * R_sh - b)
    0 < ϕ < 1 || error("Implied bargaining power ϕ = $ϕ lies outside (0, 1).")
    x_m  = Q / e^ξ_inv
    a    = repost ? Q / ((Λ_r - 1 + p_r) / p_r)^(1 / ψ_r) : 0.0

    return (δ = δ, s = s, r = r, A = A, η_L = η_L, b = b, ϕ = ϕ, κ = κ, x_m = x_m,
            ξ_inv = ξ_inv, zbar = zbar, a = a, p_r = repost ? p_r : 0.0, ψ_r = ψ_r,
            targets = tg)
end

# ── Steady state from the parameters ────────────────────────────────────────

"""
    steady_state_ags(cal) → NamedTuple

Steady state solved from the parameters alone (not from the targets): θ solves job creation;
given θ, Q solves Q = x_m e(Q)^ξ_inv with e(Q) = δ(v+1−u) + (1−δ)(1−Λ_r(Q)) s (1−u).
"""
function steady_state_ags(cal)
    @unpack δ, s, r, A, η_L, b, ϕ, κ, x_m, ξ_inv, zbar = cal
    τ = δ + (1 - δ) * s

    function given_θ(θ)
        q = A * θ^(-η_L)
        f = θ * q
        u = τ / (τ + (1 - δ) * f)
        v = θ * u
        e_of(Q) = δ * (v + 1 - u) + (1 - δ) * (1 - Λr_ags(Q, cal)) * s * (1 - u)
        # Q − x_m e(Q)^ξ_inv is increasing in Q (e falls as reactivation rises): one root
        Q = ξ_inv == 0 ? x_m : find_zero(Q -> Q - x_m * e_of(Q)^ξ_inv, (1e-12, 1e6), Bisection())
        Λ_r  = Λr_ags(Q, cal)
        Qrep = Qrep_ags(Q, Λ_r, cal)
        K    = Q * (r + δ) / (1 + r)
        R_sh = (1 - δ) * (Q - Qrep) / (1 + r)
        w    = ϕ * (zbar - K + θ * (K + q * κ) - s * R_sh) + (1 - ϕ) * b
        jcc  = zbar - w - K - s * R_sh - (r + τ) / (1 - δ) * (κ + K / q)
        return jcc, (; θ, q, f, u, v, e = e_of(Q), Q, Λ_r, Qrep, K, R_sh, w,
                       H = (1 - δ) * Qrep / (1 + r), v_pre = v - e_of(Q))
    end

    θ = find_zero(t -> given_θ(t)[1], (1e-4, 1e2), Bisection())
    return given_θ(θ)[2]
end

"""Upper bound of F_r: Λ_r = 1 once Q/Q̄ exceeds this factor (Notes/ags_equilibrium.md §5)."""
bound_factor(cal, Λ̄_r) = cal.a == 0 ? Inf : (cal.p_r / (Λ̄_r - 1 + cal.p_r))^(1 / cal.ψ_r)
