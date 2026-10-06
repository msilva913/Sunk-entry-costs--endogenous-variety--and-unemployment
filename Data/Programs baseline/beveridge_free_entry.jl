# beveridge_free_entry.jl
# =============================================================================
# The Beveridge-curve response to a destruction shock, and how it depends on
# δ_e, in two parts. Companion to Notes/predetermined_vacancies_and_the_CK_special_case.md.
#
# PART 1. The vacancy accounting identity, in consistent units.
#   On the date destruction hits, with lagged variables at steady state,
#       ṽ_t = ẽ_t − X̄ · δ̃_e,t ,     X̄ ≡ (1−q)v + s(1−u) = (v − e)/(1 − δ_e)
#   (tilde = level deviation). Vacancies rise iff ẽ_t > X̄ δ̃_e,t. The threshold
#   on the entry ELASTICITY (ê per δ̂_e, both log deviations) is therefore
#       T_log = δ_e X̄ / e = X̄/(v + L)                     [e = δ_e(v + L)]
#   and the threshold on ê per unit LEVEL shock is T_lev = X̄/e. The Oct 5 note
#   reported T_lev and read it as T_log. They differ by the factor δ_e.
#
# PART 2. Free entry (ξ → ∞) in a nested durable-vacancy DMP economy.
#   At free entry the identity cannot sign ṽ, because e is a residual: Q = x_m
#   pins the vacancy value, the JCC pins θ, and v = θu. The sign of the
#   Beveridge response is the sign of θ̂ + û. The nested economy keeps the
#   draft's timing, laws of motion and JCC (eq:jcc_eq, eq:Kdef, eq:u_lom,
#   eq:v_lom) and shuts four things off:
#     risk neutrality (m = β), no variety channel (w_int fixed),
#     exogenous exit (Λ ≡ 1, so δ_e,t+1 = δ_t), costless reposting (Λ_r = 1).
#   Each δ_e is calibrated exactly as calibrate_shares PATH B does it, with the
#   mechanism-runner targets (b_ratio = 0.9, x_v = 0.5, dest_elast_target = 5).
#   The free-entry economy shares that steady state: set x_m = Q̄, then
#   K̄ = Q̄(r+δ_e)/(1+r) whatever ξ is.
#
# TIMING (draft eq:delta_e_lom, eq:Kdef). The innovation to δ_t at date 0 raises
#   δ_e,1. At date 0 only prices move: the JCC prices next period's destruction.
#   Positions and matches are destroyed at date 1.
#
#   Let d_t ≡ δ̃_e,t+1 (the destruction rate known at t), d_t = ρ^t d_0.
#   Guess θ̂_t = a·d_t. With Q pinned, K_t = Q[1 − β(1 − δ_e,t+1)], so
#       K̂_t = d_t/(r + δ_e)                                     (user cost)
#   The JCC C_t = β(1−δ_e,t+1) S_{t+1}, C ≡ κ + K/q, S ≡ (1−ϕ)(w_int−K−b) + (1−s−ϕf)C
#   log-linearizes to a·Den = −(F1 + F2 + F3), with B ≡ β(1−δ_e), χ_K ≡ x_v:
#       F1 = 1/(1−δ_e)                                  match survival
#       F2 = χ_K/(r+δ_e) · [1 − Bρ(1−s−ϕf)]             user cost of the slot,
#                                                        net of future hiring cost saved
#       F3 = Bρ(1−ϕ)(K̄/C̄)/(r+δ_e)                       user cost eroding match surplus
#       Den = χ_K η + Bρ[ϕ f(1−η) − (1−s−ϕf) χ_K η]
#   A flow-cost DMP (K fixed, as in standard free entry) keeps F1 only.
#
# SELF-CONTAINED (Base + Printf). The Julia packages steady_state.jl needs are not
# installed on this machine, so Stages 1–5 of calibrate_shares PATH B are mirrored
# below and checked by assertion (the steady-state JCC must hold).
#
# OUTPUT  beveridge_free_entry.txt
# =============================================================================

using Printf

# --- Targets, mirrored from steady_state.jl TARGETS + mechanism runners ------
const TG = (X_Y = 0.015, f = 0.41, η_L = 0.6, q = 0.80, sep = 0.031, r_ann = 0.04,
            ε = 4.3, N = 1.0, w = 1.0, dest_end_frac = 0.5, p_0 = 0.5,
            b_ratio = 0.9, x_v = 0.5, ξ_inv = 1.0, dest_elast_target = 5.0)
const ρ_δ = 0.592      # run_solution.jl:39, monthly
const σ_δ = 0.0669     # run_solution.jl:40, log innovation to δ

δ_e_from_ann(a) = 1 - (1 - a)^(1 / 12)

const CALIBS = [
    ("D1 settled (3.20%/yr)",     δ_e_from_ann(0.0320)),
    ("Code baseline (7.54%/yr)",  δ_e_from_ann(0.0754)),
    ("BGM/GS (10%/yr)",           δ_e_from_ann(0.10)),
    ("CK-like (δ_e = τ)",         TG.sep * (1 - 1e-9)),   # s → 0; ζ_ss stays interior
]

"Bisection on a sign change; scans a log grid for the bracket."
function root(g, lo, hi; n = 400)
    xs = exp.(range(log(lo), log(hi), length = n))
    i = findfirst(k -> sign(g(xs[k])) != sign(g(xs[k+1])), 1:n-1)
    isnothing(i) && error("no sign change in [$lo, $hi]")
    a, b = xs[i], xs[i+1]
    for _ in 1:200
        m = (a + b) / 2
        sign(g(m)) == sign(g(a)) ? (a = m) : (b = m)
    end
    return (a + b) / 2
end

# --- calibrate_shares, PATH B, Stages 1–5 (steady_state.jl:669–840) ---------
function calibrate(δ_e; t = TG)
    δ   = (1 - t.dest_end_frac) * δ_e
    surv = (1 - δ_e) / (1 - δ)
    μ   = t.ε / (t.ε - 1)
    τ   = t.sep
    s   = (τ - δ_e) / (1 - δ_e)
    r   = (1 + t.r_ann)^(1 / 12) - 1
    f   = t.f / (1 - δ_e);  q = t.q / (1 - δ_e)        # conditional-on-survival rates
    θ   = f / q
    u   = τ / (τ + (1 - δ_e) * f)
    v   = θ * u;  L = 1 - u
    e   = δ_e * (v + L)
    N_e = δ_e / (1 - δ_e) * t.N
    b   = t.b_ratio * t.w
    ρ   = t.N^(1 / (t.ε - 1))
    ζ   = (surv - (1 - t.p_0)) / t.p_0
    ψ   = t.dest_elast_target * (1 - ζ) / ζ
    cons = ψ / (1 + ψ) * t.p_0
    π_s = (μ - 1) / μ * (1 - cons) * (r + δ_e) / (r + δ_e + cons * (1 - δ_e))
    L_c = (r + δ_e) * L / (r + δ_e + δ_e * π_s * μ)
    sr  = (r + τ) / (1 - δ_e) / (q * t.x_v)

    function stage4(Q)
        K     = Q * (r + δ_e) / (1 + r)
        κ     = (1 - t.x_v) / t.x_v * K / q
        X     = e / (1 + t.ξ_inv) * Q + κ * q * v
        w_int = sr * K + t.w + K
        z     = (μ / ρ) * w_int
        f_e   = π_s * z * L_c * (1 - δ_e) * μ / (t.N * (r + δ_e))
        ν_f   = ρ * f_e / μ
        Y_c   = ρ * z * L_c
        x_c   = Y_c / (t.ε * t.N) + ν_f
        C     = Y_c - X - t.N * cons * x_c
        Y     = C + ν_f * N_e
        return 100 * (X / t.X_Y - Y) / (X / t.X_Y + Y), (; K, κ, w_int)
    end
    Q = root(x -> stage4(x)[1], 1e-3, 1e5)
    _, o = stage4(Q)
    ϕ = (t.w - b) / (o.w_int - o.K + θ * (o.K + q * o.κ) - b)

    # Steady-state JCC (steady_state.jl:358) must hold at the calibrated point
    lhs = (o.κ + o.K / q) * (r + τ + (1 - δ_e) * ϕ * q * θ)
    rhs = (1 - δ_e) * (1 - ϕ) * (o.w_int - o.K - b)
    @assert isapprox(lhs, rhs; rtol = 1e-8) "JCC fails at δ_e = $δ_e"

    return (; δ_e, τ, s, r, f, q, θ, u, v, L, e, b, ϕ, Q, o.K, o.κ, o.w_int,
              C = o.κ + o.K / q, χ_K = (o.K / q) / (o.κ + o.K / q), η = t.η_L)
end

# --- Part 2: free-entry tightness coefficient a, with θ̂_t = a·d_t ----------
function tightness(c; ρ = ρ_δ, user_cost = true)
    B   = (1 - c.δ_e) / (1 + c.r)
    m   = 1 - c.s - c.ϕ * c.f                       # continuation weight on C_{t+1}
    F1  = 1 / (1 - c.δ_e)
    F2  = user_cost ? c.χ_K / (c.r + c.δ_e) * (1 - B * ρ * m) : 0.0
    F3  = user_cost ? B * ρ * (1 - c.ϕ) * (c.K / c.C) / (c.r + c.δ_e) : 0.0
    Den = c.χ_K * c.η + B * ρ * (c.ϕ * c.f * (1 - c.η) - m * c.χ_K * c.η)
    return (; a = -(F1 + F2 + F3) / Den, F1, F2, F3, Den)
end

# Steady-state (ρ = 1) comparative static, solved directly from the SS JCC,
# as a check on the ρ = 1 case of `tightness`.
function tightness_permanent(c)
    Ω  = c.r + c.τ + (1 - c.δ_e) * c.ϕ * c.f
    dK = c.K / (c.r + c.δ_e)                         # dK/dδ_e with Q pinned
    # d/dδ_e of  C·Ω − (1−δ_e)(1−ϕ)(w_int−K−b) = 0, with C = κ + K/q(θ)
    dθ_part = c.χ_K * c.C * c.η * Ω + c.C * (1 - c.δ_e) * c.ϕ * c.f * (1 - c.η)
    dδ_part = dK / c.q * Ω + c.C * (1 - c.s - c.ϕ * c.f) +
              (1 - c.ϕ) * (c.w_int - c.K - c.b) + (1 - c.δ_e) * (1 - c.ϕ) * dK
    return -dδ_part / dθ_part
end

# --- Part 2: impulse responses of u, v, e under free entry ------------------
# Shock: d_0 = 1 (one unit of monthly δ_e, level), δ̃_e,t = ρ^(t−1) for t ≥ 1.
function irf(c, a; ρ = ρ_δ, H = 48)
    Dt(t) = t >= 1 ? ρ^(t - 1) : 0.0               # δ̃_e,t
    θh = [a * ρ^t for t in 0:H]                     # θ̂_t
    ũ = zeros(H + 1); ṽ = zeros(H + 1); ẽ = zeros(H + 1)
    X̄ = (1 - c.q) * c.v + c.s * (1 - c.u)           # bracket in eq:v_lom
    for t in 0:H
        i = t + 1
        if t == 0
            ũ[i] = 0.0                              # u_0 predetermined
        else
            ũ[i] = (1 - (1 - c.δ_e) * c.f - c.τ) * ũ[i-1] +
                   (c.u * c.f + (1 - c.u) * (1 - c.s)) * Dt(t) -
                   c.u * (1 - c.δ_e) * c.f * (1 - c.η) * θh[i-1]
        end
        ṽ[i] = c.v * θh[i] + c.θ * ũ[i]             # v = θu
        ṽpre = t == 0 ? 0.0 :
               -X̄ * Dt(t) + (1 - c.δ_e) * ((1 - c.q) * ṽ[i-1] +
               c.v * c.η * c.q * θh[i-1] - c.s * ũ[i-1])
        ẽ[i] = ṽ[i] - ṽpre                          # entry is the residual
    end
    return (; ũ, ṽ, ẽ)
end


# =============================================================================
lines = String[]
P(s) = push!(lines, s)
P("="^100)
P("THE BEVERIDGE RESPONSE TO A DESTRUCTION SHOCK: ACCOUNTING, AND FREE ENTRY")
P("Targets: f=0.41 q=0.80 sep=0.031 X_Y=0.015 b_ratio=0.9 x_v=0.5 dest_elast_target=5 (PATH B)")
P(@sprintf("Shock: ρ_δ = %.3f monthly; 1-SD innovation σ_δ = %.4f (log δ) ≈ %.1f%% of δ̄_e in levels",
           ρ_δ, σ_δ, 100 * (1 - TG.dest_end_frac) * σ_δ))
P("="^100)
P("")
P("PART 1. Vacancy accounting on the destruction date (lagged variables at steady state)")
P("  ṽ = ẽ − X̄·δ̃_e.   T_log: needed ê per 1% rise in δ_e.   T_lev: needed ê per unit level δ̃_e.")
P("  X̄/v: share of the stock exposed.   hole: vacancies destroyed by a 1% rise in δ_e, % of v̄.")
P("")
P(@sprintf("%-28s %9s %8s %8s %8s %10s %10s %10s", "Calibration", "δ_e", "δ_e/τ", "e/v",
           "X̄/v", "T_log", "T_lev", "hole"))
P("-"^100)
cals = [(n, calibrate(d)) for (n, d) in CALIBS]
for (n, c) in cals
    X̄ = (1 - c.q) * c.v + c.s * (1 - c.u)
    @assert isapprox(X̄ * (1 - c.δ_e), c.v - c.e; rtol = 1e-10)   # SS vacancy LOM
    T_log = c.δ_e * X̄ / c.e
    @assert isapprox(T_log, X̄ / (c.v + c.L); rtol = 1e-10)
    hole = 100 * 0.01 * c.δ_e * X̄ / c.v
    P(@sprintf("%-28s %9.6f %8.4f %8.4f %8.4f %10.4f %10.2f %9.4f%%", n, c.δ_e, c.δ_e / c.τ,
               c.e / c.v, X̄ / c.v, T_log, X̄ / c.e, hole))
end
P("-"^100)
P("")
P("PART 1b. Effective within-period elasticity of the vacancy STOCK to Q: ∂v̂/∂Q̂ = ξ·e/v")
P("  (entry is the only within-period margin under costless reposting). CK estimate ξ = 0.265.")
P("")
P(@sprintf("%-28s %10s %10s %10s %10s %14s", "Calibration", "ξ = 0.265", "ξ = 1", "ξ = 2",
           "ξ = 10", "ξ for unit"))
P("-"^100)
for (n, c) in cals
    ev = c.e / c.v
    P(@sprintf("%-28s %10.3f %10.3f %10.3f %10.3f %14.1f", n, 0.265ev, ev, 2ev, 10ev, 1 / ev))
end
P("-"^100)
P("")
P("PART 2. Free entry (Q = x_m) in the nested durable-vacancy economy")
P("  a = θ̂ per unit level of δ̃_e,t+1. Forces: F1 match survival, F2 user cost of the slot,")
P("  F3 user cost eroding the match surplus. 'DMP' = flow-cost vacancy (K fixed): F1 only.")
P("")
P(@sprintf("%-28s %7s %7s %7s %8s %8s %8s %8s %9s %9s %9s", "Calibration", "ϕ", "χ_K", "K/C",
           "F1", "F2", "F3", "Den", "a", "a_DMP", "(r+δ_e)/τ"))
P("-"^100)
for (n, c) in cals
    tz = tightness(c)
    a_dmp = tightness(c; user_cost = false).a
    @assert isapprox(tightness(c; ρ = 1.0).a, tightness_permanent(c); rtol = 1e-6) "ρ=1 check fails: $n"
    P(@sprintf("%-28s %7.4f %7.3f %7.3f %8.2f %8.2f %8.2f %8.3f %9.1f %9.2f %9.3f", n, c.ϕ, c.χ_K,
               c.K / c.C, tz.F1, tz.F2, tz.F3, tz.Den, tz.a, a_dmp, (c.r + c.δ_e) / c.τ))
end
P("-"^100)
P("")
P("  Impulse responses at free entry, per 1% rise in δ_e,1 (d_0 = 0.01·δ̄_e). Levels, % of SS.")
P("  t=0: anticipation (no destruction yet). t=1: destruction date. slope = Σṽũ/Σũ² over 48 months.")
P("  e≥0 limit: largest shock, in % of δ̄_e, before entry would have to turn negative.")
P("  v̂>0 @: first month in which vacancies are above steady state.")
P("")
P(@sprintf("%-28s %-6s %9s %9s %9s %9s %9s %10s %7s", "Calibration", "model", "v̂_0", "û_1",
           "v̂_1", "û_peak", "slope", "e≥0 limit", "v̂>0 @"))
P("-"^100)
for (n, c) in cals, (lab, uc) in (("dur", true), ("DMP", false))
    a = tightness(c; user_cost = uc).a
    R = irf(c, a)
    sc = 0.01 * c.δ_e                                   # 1% of δ̄_e
    v̂ = 100 .* R.ṽ .* sc ./ c.v
    û = 100 .* R.ũ .* sc ./ c.u
    slope = sum(R.ṽ .* R.ũ) / sum(R.ũ .^ 2)
    emin = minimum(R.ẽ)
    lim = emin < 0 ? 100 * (c.e / -emin) / c.δ_e : Inf
    tpos = findfirst(>(0), v̂)
    P(@sprintf("%-28s %-6s %8.4f%% %8.4f%% %8.4f%% %8.4f%% %9.4f %9.1f%% %7s", n, lab, v̂[1], û[2],
               v̂[2], maximum(û), slope, lim, isnothing(tpos) ? "never" : string(tpos - 1)))
end
P("-"^100)
P("")
P("READING")
P("  dur = durable vacancy at free entry (this model, ξ → ∞). DMP = same SS, flow-cost vacancy.")
P("  slope < 0 is a downward-sloping Beveridge response.")
P("")

# --- Part 3: the destruction-date sign condition and its robustness --------
# On date 1, θ̂_1 = aρd_0 and û_1 = [c_u − ū(1−δ_e)f(1−η)a]d_0/ū, c_u ≡ ūf + (1−ū)(1−s), so
#     v̂_1/d_0 = a·[ρ − (1−δ_e)f(1−η)] + c_u/ū .
# For ρ > (1−δ_e)f(1−η), vacancies fall on the destruction date iff |a| > a* ≡ (c_u/ū)/[ρ − (1−δ_e)f(1−η)].
function sign_condition(c, a; ρ = ρ_δ)
    c_u   = c.u * c.f + (1 - c.u) * (1 - c.s)
    g     = ρ - (1 - c.δ_e) * c.f * (1 - c.η)
    v1    = a * g + c_u / c.u
    R = irf(c, a; ρ = ρ)
    @assert isapprox(v1 * c.v, R.ṽ[2]; rtol = 1e-8)        # closed form = simulated IRF
    return (; astar = (c_u / c.u) / g, v1, g,
              slope = sum(R.ṽ .* R.ũ) / sum(R.ũ .^ 2))
end

P("PART 3. Free entry: destruction-date sign condition  v̂_1/d_0 = a·[ρ − (1−δ_e)f(1−η)] + c_u/ū")
P("  Vacancies fall on date 1 iff |a| > a*. 'slope' is the 48-month path slope as above.")
P("")
P(@sprintf("%-28s %9s %9s %9s %10s %9s", "Calibration", "|a|", "a*", "|a|/a*", "v̂_1/d_0", "slope"))
P("-"^100)
for (n, c) in cals
    a = tightness(c).a
    sc = sign_condition(c, a)
    P(@sprintf("%-28s %9.1f %9.1f %9.2f %10.1f %9.4f", n, -a, sc.astar, -a / sc.astar,
               sc.v1, sc.slope))
end
P("-"^100)
P("")
P("  Robustness at free entry: D1 and CK-like, varying one input at a time from the")
P("  mechanism targets (x_v = 0.5, b_ratio = 0.9, ρ_δ = 0.592). x_v = sunk (durable) share")
P("  of recruiting cost; x_v → 0 is a pure flow-cost vacancy.")
P("")
P(@sprintf("%-14s %-8s %9s %9s %9s %10s %10s", "varied", "value", "|a| D1", "a* D1",
           "slope D1", "|a| CK", "slope CK"))
P("-"^100)
δD1 = δ_e_from_ann(0.0320); δCK = TG.sep * (1 - 1e-9)
rows = vcat([("x_v", x, (TG..., x_v = x), ρ_δ) for x in (0.05, 0.1, 0.25, 0.5, 0.75, 1.0)],
            [("b_ratio", b, (TG..., b_ratio = b), ρ_δ) for b in (0.71, 0.8, 0.9, 0.95)],
            [("ρ_δ", ρ, TG, ρ) for ρ in (0.3, 0.592, 0.9, 0.99)])
for (lab, val, tg, ρ) in rows
    out = map((δD1, δCK)) do δ
        try
            c = calibrate(δ; t = tg)
            a = tightness(c; ρ = ρ).a
            (a, sign_condition(c, a; ρ = ρ))
        catch
            nothing
        end
    end
    f(x, k) = isnothing(x) ? "   fail" : @sprintf("%9.4g", k(x))
    P(@sprintf("%-14s %-8.3g %s %s %s %10s %10s", lab, val,
               f(out[1], x -> -x[1]), f(out[1], x -> x[2].astar), f(out[1], x -> x[2].slope),
               f(out[2], x -> -x[1]), f(out[2], x -> x[2].slope)))
end
P("-"^100)
P("")
# Where the path slope changes sign, by bisection on one input at a time
slope_at(δ, ρ) = (c = calibrate(δ); sign_condition(c, tightness(c; ρ = ρ).a; ρ = ρ).slope)
function flip(g, lo, hi)
    sign(g(lo)) == sign(g(hi)) && return NaN
    for _ in 1:60
        m = (lo + hi) / 2
        sign(g(m)) == sign(g(lo)) ? (lo = m) : (hi = m)
    end
    return (lo + hi) / 2
end
P("  Sign flips of the path slope at free entry:")
for (n, d) in CALIBS[1:3]
    ρc = flip(ρ -> slope_at(d, ρ), 0.05, 0.99)
    P(@sprintf("    %-28s slope < 0 iff ρ_δ > %.3f (monthly; %.3f quarterly)", n, ρc, ρc^3))
end
ann(δ) = 1 - (1 - δ)^12
δc = flip(δ -> slope_at(δ, ρ_δ), δ_e_from_ann(0.0320), TG.sep * (1 - 1e-9))
P(@sprintf("    at ρ_δ = %.3f: slope < 0 iff δ_e < %.5f monthly (%.1f%%/yr, δ_e/τ = %.3f)",
           ρ_δ, δc, 100 * ann(δc), δc / TG.sep))
ρc_hi = flip(ρ -> slope_at(TG.sep * (1 - 1e-9), ρ), 0.05, 0.999)
P(isnan(ρc_hi) ? "    CK-like: slope > 0 for every ρ_δ in [0.05, 0.999]" :
  @sprintf("    CK-like: slope < 0 iff ρ_δ > %.3f", ρc_hi))
P("="^100)

text = join(lines, "\n")
println(text)
open(joinpath(@__DIR__, "beveridge_free_entry.txt"), "w") do io
    write(io, text * "\n")
end
