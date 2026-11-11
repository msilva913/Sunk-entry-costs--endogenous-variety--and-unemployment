# ags_steady_state.jl
# =============================================================================
# Calibration and steady state of the AGS model with costly reactivation (Part 1 of
# decision D15; equations in Notes/ags_equilibrium.md). Standalone: no variety, no product-line
# entry block, no endogenous exit. Destroyed product lines are replaced by clones at no cost.
#
# Economy (monthly). A representative household with CRRA utility C^(1−σ)/(1−σ) owns the firms
#   (shares in a mutual fund, fixed supply 1) and can trade a riskless bond in zero net supply;
#   firms value jobs and positions with its discount factor m_{t+1} = β(C_{t+1}/C_t)^(−σ). In
#   steady state m = β = 1/(1+r) whatever σ, so the calibration and the steady state do not
#   depend on σ; σ = 0 is linear utility (Gabrovski–Silva). See Notes/ags_equilibrium.md §2.
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
Reactivation: p_r = 0 with Λr_target = 1 is costless reactivation (Λ_r ≡ 1, Q^rep = Q:
Gabrovski–Silva with three shocks); otherwise 1 − p_r < Λr_target < 1. The defaults are the D15
working point (Λ̄_r = 0.25, p_r = 0.8, ψ_r = 1).
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
    Λr_target = 0.25,     # steady-state reactivation rate (1 when p_r = 0)
    p_r       = 0.8,      # mass of the costly part of F_r; 0 = costless reactivation
    ψ_r       = 1.0,      # shape of F_r's power law
    σ         = 1.0,      # inverse IES (CRRA; log utility, as the full model's TARGETS); 0 =
                          # linear utility (the Gabrovski–Silva benchmark). Dynamics only.
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

# Reactivation-cost distribution F_r (Notes/ags_equilibrium.md §2): c = 0 with probability
# 1 − p_r; with probability p_r, c ~ power law on [0, a]: Pr(c ≤ x) = (x/a)^ψ_r. A vacated
# position is reactivated iff c ≤ Q. Both functions below are exact for every Q > 0; first-order
# perturbation (f[9], f[10] in ags_run_solution.jl) uses the interior branch Q ≤ a.

"""
    Λr_ags(Q, cal)
Reactivation rate Λ_r = F_r(Q): (1 − p_r) + p_r (Q/a)^ψ_r while Q ≤ a; 1 once Q > a (every
costly draw is accepted). Identically 1 when p_r = 0.
"""
function Λr_ags(Q, cal)
    @unpack p_r, ψ_r, a = cal
    if Q <= a
        return (1 - p_r) + p_r * (Q / a)^ψ_r
    else
        return 1.0
    end
end

"""
    Qrep_ags(Q, cal)
Expected net value of a vacated position, E max(Q − c, 0) = Q·Λ_r − E[c; c ≤ Q]. With
ψ_c,r = ψ_r/(ψ_r + 1), the mean of a power-law draw on [0, x] relative to x:
- Q ≤ a: accepted costly draws (share Λ_r − 1 + p_r) cost ψ_c,r·Q on average, so
  Q^rep = Q[Λ_r − ψ_c,r (Λ_r − 1 + p_r)];
- Q > a: every draw is accepted (Λ_r = 1) and the costly ones cost ψ_c,r·a on average, so
  Q^rep = Q − p_r ψ_c,r a.
Q when p_r = 0.
"""
function Qrep_ags(Q, cal)
    @unpack p_r, ψ_r, a = cal
    ψ_cr = ψ_r / (ψ_r + 1)
    if Q <= a
        Λ_r = Λr_ags(Q, cal)
        return Q * (Λ_r - ψ_cr * (Λ_r - 1 + p_r))
    else
        return Q - p_r * ψ_cr * a
    end
end

# ── Calibration ─────────────────────────────────────────────────────────────

"""
    calibrate_ags(tg) → NamedTuple

Parameters of the AGS model from targets `tg` (see AGS_TARGETS). Closed form; no root-finding.
"""
function calibrate_ags(tg = AGS_TARGETS)
    @unpack dest_ann, sep, f, q, η_L, r_ann, w, b_ratio, x_v, X_Y, ξ_inv, Λr_target, p_r, ψ_r, σ = tg
    # p_r = 0: costless reactivation, so Λ̄_r must be 1. p_r > 0: Λ̄_r strictly inside
    # (1 − p_r, 1): above the free atom's floor and off the kink of min{·, 1}.
    if p_r == 0
        if Λr_target != 1
            error("p_r = 0 is costless reactivation: set Λr_target = 1 (got $Λr_target).")
        end
    else
        if !(0 < p_r ≤ 1)
            error("p_r = $p_r must lie in [0, 1].")
        end
        if !(1 - p_r < Λr_target < 1)
            error("Λr_target = $Λr_target must lie inside (1−p_r, 1) = ($(1-p_r), 1).")
        end
    end
    Λ_r = Λr_target

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

    # Steady-state ratios per unit of Q (suffix _Q): K_Q = K/Q (user cost), κ_Q = κ/Q (per-hire
    # cost, fixed relative to K/q by x_v), R_Q = 𝓡/Q (reactivation shortfall). They let the
    # X/Y target pin Q in closed form.
    K_Q  = (r + δ) / (1 + r)
    κ_Q  = (1 - x_v) / x_v * K_Q / qc
    gap  = (1 - Λ_r) + ψ_r / (ψ_r + 1) * (Λ_r - 1 + p_r)     # 0 when p_r = 0, Λ_r = 1
    R_Q  = (1 - δ) * gap / (1 + r)
    # Job creation in steady state: z̄ − w = K + s𝓡 + (r+τ)/(1−δ)·(κ + K/q) = zw_Q·Q
    zw_Q = K_Q + s * R_Q + (r + τ) / (1 - δ) * (κ_Q + K_Q / qc)
    # Recruiting costs X = eQ/(1+ξ_inv) + κ q v = X_Q·Q
    X_Q  = e / (1 + ξ_inv) + κ_Q * qc * v
    # X/Y target with Y = z̄L = (w + zw_Q·Q)L:  X_Q·Q = X_Y (w + zw_Q·Q) L
    den  = X_Q - X_Y * zw_Q * L
    if den <= 0
        error("X_Y = $X_Y is infeasible: recruiting costs cannot reach that share of output.")
    end
    Q    = X_Y * w * L / den

    K    = K_Q * Q
    κ    = κ_Q * Q
    R_sh = R_Q * Q
    zbar = w + zw_Q * Q
    b    = b_ratio * w
    ϕ    = (w - b) / (zbar - K + θ * (K + qc * κ) - s * R_sh - b)
    if !(0 < ϕ < 1)
        error("Implied bargaining power ϕ = $ϕ lies outside (0, 1).")
    end
    x_m  = Q / e^ξ_inv
    # Scale of F_r from Λ̄_r. With p_r = 0, a has no effect; Q is an inert placeholder.
    if p_r == 0
        a = Q
    else
        a = Q / ((Λ_r - 1 + p_r) / p_r)^(1 / ψ_r)
    end

    return (δ = δ, s = s, r = r, A = A, η_L = η_L, b = b, ϕ = ϕ, κ = κ, x_m = x_m,
            ξ_inv = ξ_inv, zbar = zbar, a = a, p_r = p_r, ψ_r = ψ_r, σ = σ,
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
        if ξ_inv == 0
            Q = x_m
        else
            Q = find_zero(Q -> Q - x_m * e_of(Q)^ξ_inv, (1e-12, 1e6), Bisection())
        end
        Λ_r  = Λr_ags(Q, cal)
        Qrep = Qrep_ags(Q, cal)
        K    = Q * (r + δ) / (1 + r)
        R_sh = (1 - δ) * (Q - Qrep) / (1 + r)
        w    = ϕ * (zbar - K + θ * (K + q * κ) - s * R_sh) + (1 - ϕ) * b
        jcc  = zbar - w - K - s * R_sh - (r + τ) / (1 - δ) * (κ + K / q)
        return jcc, (; θ, q, f, u, v, e = e_of(Q), Q, Λ_r, Qrep, K, R_sh, w,
                       H = (1 - δ) * Qrep / (1 + r), v_pre = v - e_of(Q))
    end

    θ = find_zero(t -> given_θ(t)[1], (1e-4, 1e2), Bisection())
    st = given_θ(θ)[2]

    # Household side. Output z̄(1 − u); investment = creating positions and hiring (X) plus
    # reactivation costs (X^r, paid on the Ṽ^s = (1−δ)s(1−u) vacated positions). b is a
    # non-market flow: it enters bargaining, not resources (as in the full model).
    Vs   = (1 - δ) * s * (1 - st.u)
    X    = st.e * st.Q / (1 + ξ_inv) + κ * st.q * st.v
    X_r  = Vs * (st.Q * st.Λ_r - st.Qrep)
    C    = zbar * (1 - st.u) - X - X_r                  # resource constraint
    D    = (zbar - st.w) * (1 - st.u) - X - X_r          # dividends: profits net of investment
    if C <= 0 || D <= 0
        error("steady_state_ags: consumption $C or dividends $D not positive.")
    end
    V    = D / r                                         # share price: V = β(V + D)
    return (; st..., X, X_r, C, D, V, Rf = 1 + r)        # Rf: gross riskless rate, 1/β
end

"""Upper bound of F_r: Λ_r = 1 once Q/Q̄ exceeds this factor (Notes/ags_equilibrium.md §5)."""
function bound_factor(cal, Λ̄_r)
    if cal.p_r == 0
        return Inf
    else
        return (cal.p_r / (Λ̄_r - 1 + cal.p_r))^(1 / cal.ψ_r)
    end
end
