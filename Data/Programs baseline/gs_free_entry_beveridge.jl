# gs_free_entry_beveridge.jl
# =============================================================================
# δ shock at free entry (ξ → ∞) in the AGS economy of prop:ags with σ = κ = 0, which
# rem:nesting identifies as Gabrovski-Silva with a time-varying δ_t: the labor-market
# block closes on its own. Draft timing (S10): δ_t is known at the start of t and acts
# at the end of t, so u_t and v_pre,t are predetermined and u_{t+1}, v_pre,t+1 carry
# the destruction. Numbers for Notes/beveridge_free_entry_GS.tex (task R15).
# NOTE: internal formulas below use LEVEL shocks d_t; the note uses the log shock
# δ_t = δ̄ exp(δ̃_t), whose coefficient is δ̄ × the level one (a_note() checks this).
#
# Log-linear system (d_t ≡ δ_t − δ̄, d_t = ρ^t d_0; hats are log deviations).
#   Free entry, eq:Kdef:   K_t = x_m[1 − β(1 − δ_t)]   ⇒   K̂_t = d_t/(r + δ̄)
#   eq:jcc_eq:  C_t = β(1−δ_t) E_t S_{t+1},  C ≡ κ + K/q,
#               S ≡ (1−ϕ)(z/μ − b − K) + (1 − s − ϕf)C,   χ ≡ (K̄/q̄)/C̄ = x_v
#   ⇒ θ̂_t = a d_t,  a = −{1/(1−δ̄) + χ[Ψ + Bρ(1−ϕ)q̄]/(r+δ̄)} / [χηΨ + Bρϕf̄(1−η)],
#     B ≡ β(1−δ̄),  Ψ ≡ 1 − Bρ(1 − s̄ − ϕf̄),  unique iff |Φ| < 1, Φ ≡ B[m − ϕf̄(1−η)/(χη)].
#   eq:u_lom:   ũ_{t+1} = λũ_t + c_u d_t − ū ω θ̂_t,
#               λ ≡ 1 − (1−δ̄)f̄ − τ̄,  c_u ≡ (1−s̄)(1−ū) + ūf̄,  ω ≡ (1−δ̄)f̄(1−η)
#   v_t = θ_t u_t (entry is the residual).
# Results checked here: h = 0 and h = 1 responses, the h = 1 sign condition
#   |a|(ρ − ω) > f̄/τ̄   (uses c_u/ū = f̄/τ̄, an identity at the steady state),
# and the path slope Σṽũ/Σũ² = θ̄ − (v̄|a|/κ_u)·ρ(1−λ²)/(1+ρλ),  κ_u ≡ c_u + ūω|a|.
#
# Parameters: our calibration (calibrate() in beveridge_free_entry.jl, mirroring
# steady_state.jl PATH B), mechanism-runner targets (b_ratio 0.9, x_v 0.5), and x_v = 1.
# The free-entry economy shares the calibrated steady state (x_m = Q̄).
#
# OUTPUT  gs_free_entry_beveridge.txt
# =============================================================================

include(joinpath(@__DIR__, "beveridge_free_entry.jl"))   # calibrate, tightness, irf, flip
using Printf

"Closed forms of the log-linear free-entry system. Kfixed = flow-cost comparator."
function gs(c; ρ = ρ_δ, Kfixed = false)
    η, δ, s, ϕ, f, q, u, τ, χ = c.η, c.δ_e, c.s, c.ϕ, c.f, c.q, c.u, c.τ, c.χ_K
    B = (1 - δ) / (1 + c.r);  m = 1 - s - ϕ * f
    Ψ = 1 - B * ρ * m
    uc = Kfixed ? 0.0 : χ * (Ψ + B * ρ * (1 - ϕ) * q) / (c.r + δ)   # user cost of the vacancy
    a  = -(1 / (1 - δ) + uc) / (χ * η * Ψ + B * ρ * ϕ * f * (1 - η))
    Φ  = B * (m - ϕ * f * (1 - η) / (χ * η))                         # forward root
    λ  = 1 - (1 - δ) * f - τ
    c_u = (1 - s) * (1 - u) + u * f
    ω  = (1 - δ) * f * (1 - η)
    κu = c_u + u * ω * (-a)
    v1 = a * (ρ - ω) + f / τ                       # v̂_{t+1}/d_0
    slope = c.θ + c.v * a * ρ * (1 - λ^2) / (κu * (1 + ρ * λ))
    Dρ = χ * η * Ψ + B * ρ * ϕ * f * (1 - η)
    # Sufficient conditions in the form (r+δ)/τ < Γ, dropping the survival term in a
    Γ1 = χ * (Ψ + B * ρ * (1 - ϕ) * q) * (ρ - ω) / (f * Dρ)
    ΓP = χ * (Ψ + B * ρ * (1 - ϕ) * q) * (ρ * (1 - λ^2) - ω * (1 + ρ * λ)) / (f * Dρ * (1 + ρ * λ))
    # Note's form: R = unemployment loading / tightness loading (unit-free; same in level or
    # log shock units), Λρ = persistence discounted by unemployment's own persistence.
    R  = (f / τ) / (-a)
    Λρ = ρ * (1 - λ^2) / (1 + ρ * λ)
    return (; a, Φ, λ, c_u, ω, κu, v1, slope, uc, Γ1, ΓP, R, Λρ)
end

"""
Tightness coefficient per unit of the log shock δ̃ (δ_t = δ̄ exp(δ̃_t)), in the form of
Notes/beveridge_free_entry_GS.tex eq:a: written with 1−τ̄ and the effective matching rates
(1−δ̄)f̄, (1−δ̄)q̄. Valid for κ = 0 (x_v = 1), the case the note treats.
"""
function a_note(c; ρ = ρ_δ)
    η, δ, ϕ, f, q, τ, r = c.η, c.δ_e, c.ϕ, c.f, c.q, c.τ, c.r
    ω_ρ = 1 + r - ρ * (1 - τ)                       # r + τ̄ + (1−ρ)(1−τ̄)
    num = (1 + r) * δ / (1 - δ) +
          δ / (r + δ) * (ω_ρ + ρ * ϕ * (1 - δ) * f + ρ * (1 - ϕ) * (1 - δ) * q)
    den = η * ω_ρ + ρ * ϕ * (1 - δ) * f
    return (; a = -num / den, Φ = (1 - τ - ϕ * (1 - δ) * f / η) / (1 + r))
end

lines = String[]
P(s) = push!(lines, s)
P("="^100)
P("δ SHOCK AT FREE ENTRY UNDER prop:independence (σ = 0, ζ = 0, p_0 = 0): LOG-LINEAR DYNAMICS")
P("Draft timing (S10). d_0 = one unit of δ_t (level). v̂, û in log points per unit of d_0.")
P("="^100)

for (tag, tg) in (("x_v = 0.5 (mechanism runners)", TG), ("x_v = 1 (κ = 0)", (TG..., x_v = 1.0)))
    cs = [(n, calibrate(d; t = tg)) for (n, d) in CALIBS]
    P("")
    P("── " * tag * " " * "─"^(96 - length(tag)))
    P("")
    P(@sprintf("%-26s %6s %8s %8s %9s %9s %9s %9s %9s", "Calibration", "Φ", "a(ρ=0)",
               "a(ρ̂)", "v̂_0/d", "û_1/d", "v̂_1/d", "f/τ", "|a|(ρ−ω)"))
    P("-"^100)
    for (n, c) in cs
        g0 = gs(c; ρ = 0.0);  g = gs(c)
        @assert abs(g.Φ) < 1 "indeterminate: $n"
        @assert isapprox(g.c_u / c.u, c.f / c.τ; rtol = 1e-10)       # c_u/ū = f/τ
        @assert isapprox(g.a, tightness(c).a; rtol = 1e-10)           # same coefficient
        @assert isapprox(g0.a, tightness(c; ρ = 0.0).a; rtol = 1e-10)
        for ρ in (0.0, 0.3, ρ_δ, 0.9)                                 # proposition in terms of R
            gρ = gs(c; ρ = ρ)
            @assert isapprox(gρ.slope, c.θ * (1 - gρ.Λρ / (gρ.ω + gρ.R)); rtol = 1e-10)
            @assert (gρ.v1 < 0) == (ρ > gρ.ω + gρ.R)
            @assert (gρ.slope < 0) == (gρ.Λρ > gρ.ω + gρ.R)
        end
        if isapprox(c.χ_K, 1.0)                                       # the note's form (κ = 0)
            for ρ in (0.0, 0.3, ρ_δ, 0.9)
                @assert isapprox(a_note(c; ρ = ρ).a, c.δ_e * gs(c; ρ = ρ).a; rtol = 1e-10) "eq:a fails: $n"
                @assert isapprox(a_note(c; ρ = ρ).Φ, gs(c; ρ = ρ).Φ; rtol = 1e-10)
            end
        end
        # IRF checks (beveridge_free_entry.jl's simulation uses the same timing)
        R  = irf(c, g.a; H = 3000)
        @assert isapprox(R.ṽ[1] / c.v, g.a; rtol = 1e-10)             # h = 0: v̂ = θ̂, û = 0
        @assert R.ũ[1] == 0.0
        @assert isapprox(R.ũ[2], g.κu; rtol = 1e-10)                  # h = 1
        @assert isapprox(R.ṽ[2] / c.v, g.v1; rtol = 1e-10)
        @assert isapprox(sum(R.ṽ .* R.ũ) / sum(R.ũ .^ 2), g.slope; rtol = 1e-6)
        R0 = irf(c, g0.a; ρ = 0.0, H = 50)                            # iid: v = θ̄u from t+1
        @assert all(isapprox.(R0.ṽ[2:end], c.θ .* R0.ũ[2:end]; atol = 1e-12))
        P(@sprintf("%-26s %6.3f %8.1f %8.1f %9.1f %9.1f %9.1f %9.2f %9.1f", n, g.Φ, g0.a, g.a,
                   g.a, g.κu / c.u, g.v1, c.f / c.τ, -g.a * (ρ_δ - g.ω)))
    end
    P("-"^100)
    P("  h = 1: vacancies fall iff |a|(ρ − ω) > f/τ.  ρ = 0 (iid): v̂_1 = û_1 > 0 and v = θ̄u from t+1.")
    P("")
    P(@sprintf("%-26s %10s %10s %10s %12s %12s %10s", "Calibration", "slope iid", "slope ρ̂",
               "slope K̄", "ρ*: h=1", "ρ*: path", "(r+δ)/τ"))
    P("-"^100)
    for (n, c) in cs
        ρ1 = flip(ρ -> gs(c; ρ = ρ).v1, 1e-3, 0.999)
        ρp = flip(ρ -> gs(c; ρ = ρ).slope, 1e-3, 0.999)
        P(@sprintf("%-26s %10.4f %10.4f %10.4f %12s %12s %10.3f", n, gs(c; ρ = 0.0).slope,
                   gs(c).slope, gs(c; Kfixed = true).slope,
                   isnan(ρ1) ? "none" : @sprintf("%.3f", ρ1),
                   isnan(ρp) ? "none" : @sprintf("%.3f", ρp), (c.r + c.δ_e) / c.τ))
    end
    P("-"^100)
    P("  slope = Σṽũ/Σũ² over the full path (levels). slope K̄: flow-cost vacancy (K fixed), ρ̂.")
    P("  ρ*: persistence above which v̂_1 < 0 (h=1) or the path slope < 0 ('none' = never in (0,1)).")
    P("")
    P("  Sufficient conditions (r+δ)/τ < Γ at ρ̂ (survival term dropped from a):")
    P(@sprintf("%-26s %10s %10s %10s", "Calibration", "(r+δ)/τ", "Γ_1 (h=1)", "Γ_P (path)"))
    P("-"^100)
    for (n, c) in cs
        g = gs(c)
        P(@sprintf("%-26s %10.3f %10.3f %10.3f", n, (c.r + c.δ_e) / c.τ, g.Γ1, g.ΓP))
    end
    P("-"^100)
    P("")
    P("  First-quarter averages (months 0–2 after the shock), log points per unit of d_0:")
    P(@sprintf("%-26s %12s %12s %12s %12s", "Calibration", "ū_Q1 iid", "v̄_Q1 iid", "ū_Q1 ρ̂",
               "v̄_Q1 ρ̂"))
    P("-"^100)
    for (n, c) in cs
        q1(ρ) = (R = irf(c, gs(c; ρ = ρ).a; ρ = ρ, H = 3);
                 (sum(R.ũ[1:3]) / 3 / c.u, sum(R.ṽ[1:3]) / 3 / c.v))
        a0 = q1(0.0); a1 = q1(ρ_δ)
        P(@sprintf("%-26s %12.1f %12.1f %12.1f %12.1f", n, a0[1], a0[2], a1[1], a1[2]))
    end
    P("-"^100)
end
P("")
P("── Proposition in terms of R ≡ (f̄δ̄/τ̄)/|a|: h=1 iff ρ > ω+R; path iff ρ(1−λ²)/(1+ρλ) > ω+R ──")
P(@sprintf("%-26s %-8s %8s %8s %8s %10s %9s", "Calibration", "x_v", "ω", "R", "ω+R",
           "ρ(1−λ²)/(1+ρλ)", "slope"))
P("-"^100)
for (tag, tg) in (("0.5", TG), ("1", (TG..., x_v = 1.0))), (n, d) in CALIBS
    c = calibrate(d; t = tg); g = gs(c)
    P(@sprintf("%-26s %-8s %8.4f %8.4f %8.4f %10.4f %9.4f", n, tag, g.ω, g.R, g.ω + g.R, g.Λρ, g.slope))
end
P("-"^100)
P(@sprintf("  at ρ = %.3f", ρ_δ))
P("")
P("── Exposition tables: κ = 0 (x_v = 1), per 1% rise in δ_t (d_0 = 0.01·δ̄), in % ──────────────────")
for (n, d) in CALIBS[[1, 4]]
    c = calibrate(d; t = (TG..., x_v = 1.0))
    sc = 0.01 * c.δ_e
    P("")
    P(@sprintf("%s:  1/(r+δ̄) = %.1f   f̄/τ̄ = %.2f   λ = %.3f   ω = (1−δ̄)f̄(1−η) = %.3f   λ+ω = %.3f",
               n, 1 / (c.r + c.δ_e), c.f / c.τ, gs(c).λ, gs(c).ω, gs(c).λ + gs(c).ω))
    P(@sprintf("    ē/v̄ = %.4f   X̄/v̄ = %.4f   (1−δ̄)(1−q̄) = %.4f   (1−δ̄)ηq̄ = %.4f   (1−δ̄)s̄/θ̄ = %.4f",
               c.e / c.v, ((1 - c.q) * c.v + c.s * (1 - c.u)) / c.v, (1 - c.δ_e) * (1 - c.q),
               (1 - c.δ_e) * c.η * c.q, (1 - c.δ_e) * c.s / c.θ))
    # Log-shock form δ_t = δ̄ exp(δ̃_t): weights per unit of δ̃ are δ̄ × the level weights
    P(@sprintf("    log shock δ̃:  δ̄/(r+δ̄) = %.4f   f̄δ̄/τ̄ = %.4f   (1−s̄)δ̄/τ̄ = %.4f   (1−δ̄)f̄ = %.4f",
               c.δ_e / (c.r + c.δ_e), c.f * c.δ_e / c.τ, (1 - c.s) * c.δ_e / c.τ, (1 - c.δ_e) * c.f))
    P(@sprintf("    log shock δ̃:  a·δ̄ = %.4f (iid), %.4f (ρ = %.3f);  r = %.6f",
               c.δ_e * gs(c; ρ = 0.0).a, c.δ_e * gs(c).a, ρ_δ, c.r))
    for (lab, ρ) in (("iid", 0.0), (@sprintf("ρ = %.3f", ρ_δ), ρ_δ))
        a = gs(c; ρ = ρ).a
        R = irf(c, a; ρ = ρ, H = 12)
        P(@sprintf("    %-11s a = %.1f", lab, a))
        P(@sprintf("    %4s %9s %9s %9s %9s", "t", "θ̂_t", "û_t", "v̂_t", "ê_t"))
        for t in 0:6
            i = t + 1
            P(@sprintf("    %4d %9.3f %9.3f %9.3f %9.2f", t, 100 * a * ρ^t * sc,
                       100 * R.ũ[i] * sc / c.u, 100 * R.ṽ[i] * sc / c.v, 100 * R.ẽ[i] * sc / c.e))
        end
        tpos = findfirst(>(0), R.ṽ)
        P(@sprintf("    first month with v̂ > 0: %s", isnothing(tpos) ? "none ≤ 12" : string(tpos - 1)))
    end
end

# =============================================================================
# s shock (Oct 7, 2026). s_t = s̄ exp(s̃_t), s̃ AR(1). Free entry: K = x_m(r+δ̄)/(1+r) does
# not depend on s, so tightness moves only through the continuation term (1−s_{t+1})K/q
# of eq:jcc_eq (a match formed at t first faces separation at the end of t+1).
#   θ̂_t = b s̃_t,   b = −ρ(1−δ̄)s̄ / D,   D = η[1+r−ρ(1−τ̄)] + ρϕ(1−δ̄)f̄  (same D as a)
#   û_t = λ û_{t−1} + ℓ_s s̃_{t−1} − ω θ̂_{t−1},   ℓ_s = (1−δ̄)f̄(1−δ̄/τ̄)
#   R_s = ℓ_s/|b| = (1−δ̄)f̄ D/(ρ τ̄);  same h = 1 and path conditions as the δ shock.
# Sufficient for positive comovement at every t ≥ 1 and every ρ:
#   (1−δ̄)f̄ [η + ϕ(1−δ̄)f̄/τ̄] ≥ 1   (since R_s > that expression and both persistence terms < 1).
# Log units throughout (the note's units). κ = 0 (x_v = 1).
# =============================================================================
const ρ_s_cal = 0.874      # context/parameters.md: monthly AR(1) on s, part6b

"s-shock closed forms, per unit of the log shock s̃."
function gs_s(c; ρ = ρ_s_cal)
    η, δ, ϕ, f, τ, r, s, θ = c.η, c.δ_e, c.ϕ, c.f, c.τ, c.r, c.s, c.θ
    fe = (1 - δ) * f
    D  = η * (1 + r - ρ * (1 - τ)) + ρ * ϕ * fe
    b  = -ρ * (1 - δ) * s / D
    ℓs = fe * (1 - δ / τ)
    ω  = (1 - η) * fe
    λ  = 1 - fe - τ
    Rs = b == 0 ? Inf : ℓs / (-b)
    Λρ = ρ * (1 - λ^2) / (1 + ρ * λ)
    slope = θ * (1 - Λρ / (ω + Rs))
    bound = fe * (η + ϕ * fe / τ)
    return (; b, D, ℓs, ω, λ, Rs, Λρ, slope, bound, v1 = ℓs - (-b) * (ρ - ω))
end

"Simulate the linear s-shock system (log units, ε_0 = 1) for H periods."
function irf_s(c, g; ρ, H = 3000)
    θh = [g.b * ρ^t for t in 0:H]
    uh = zeros(H + 1)
    for t in 1:H
        uh[t+1] = g.λ * uh[t] + g.ℓs * ρ^(t - 1) - g.ω * θh[t]
    end
    return (; θh, uh, vh = θh .+ uh)
end

P("")
P("── s shock at free entry (κ = 0, log units): θ̂ = b s̃; h=1 iff ρ > ω+R_s; path iff ρ(1−λ²)/(1+ρλ) > ω+R_s ──")
P(@sprintf("%-26s %7s %9s %8s %8s %9s %9s %9s %8s", "Calibration", "ρ", "b", "ℓ_s", "R_s",
           "ω+R_s", "Λ_ρ", "slope", "bound"))
P("-"^100)
for (n, d) in CALIBS[1:3]                     # δ_e = τ has s = 0: no s margin
    c = calibrate(d; t = (TG..., x_v = 1.0))
    # Check 1: b at ρ = 1 equals the steady-state derivative d ln θ / d ln s (K fixed)
    A = c.f / c.θ^(1 - c.η); K = c.Q * (c.r + c.δ_e) / (1 + c.r)
    function θss(sv)
        τv = c.δ_e + sv * (1 - c.δ_e)
        root(th -> (K / (A * th^(-c.η))) * (c.r + τv + (1 - c.δ_e) * c.ϕ * A * th^(1 - c.η)) -
                   (1 - c.δ_e) * (1 - c.ϕ) * (c.w_int - K - c.b), 1e-3, 50.0)
    end
    h = 1e-6
    fd = (log(θss(c.s * exp(h))) - log(θss(c.s * exp(-h)))) / (2h)
    @assert isapprox(gs_s(c; ρ = 1.0).b, fd; rtol = 1e-5) "b(ρ=1) fails: $n"
    # Check 2: ℓ_s equals the one-step nonlinear response of ln u to ln s at θ = θ̄
    u1(sv) = (1 - (1 - c.δ_e) * c.f) * c.u + (c.δ_e + sv * (1 - c.δ_e)) * (1 - c.u)
    @assert isapprox(gs_s(c).ℓs, (log(u1(c.s * exp(h))) - log(u1(c.s * exp(-h)))) / (2h); rtol = 1e-5)
    # Check 3: slope formula against a simulated path; sign claims under the bound
    for ρ in (0.3, ρ_δ, ρ_s_cal, 0.99)
        g = gs_s(c; ρ = ρ); S = irf_s(c, g; ρ = ρ)
        sim = sum(c.v .* S.vh .* c.u .* S.uh) / sum((c.u .* S.uh) .^ 2)
        @assert isapprox(sim, g.slope; rtol = 1e-6) "s slope fails: $n, ρ = $ρ"
        @assert isapprox(S.vh[2], g.v1; rtol = 1e-10)
        g.bound >= 1 && @assert all(S.vh[2:end] .> 0) && g.slope > 0
        P(@sprintf("%-26s %7.3f %9.4f %8.4f %8.3f %9.3f %9.3f %9.4f %8.3f", n, ρ, g.b, g.ℓs,
                   g.Rs, g.ω + g.Rs, g.Λρ, g.slope, g.bound))
    end
    g0 = gs_s(c; ρ = 0.0); S0 = irf_s(c, g0; ρ = 0.0, H = 50)
    @assert g0.b == 0 && all(isapprox.(S0.vh, S0.uh; atol = 1e-14))   # iid: θ never moves
end
P("-"^100)
P("  bound = (1−δ̄)f̄[η + ϕ(1−δ̄)f̄/τ̄] ≥ 1 ⇒ v̂_t > 0 for all t ≥ 1 and slope > 0 for every ρ (asserted).")
P(@sprintf("  For comparison, R_δ at ρ = %.3f (δ shock, κ = 0): D1 %.4f, δ_e = τ %.4f.", ρ_δ,
           gs(calibrate(CALIBS[1][2]; t = (TG..., x_v = 1.0))).R,
           gs(calibrate(CALIBS[4][2]; t = (TG..., x_v = 1.0))).R))
P("")
cD1 = calibrate(CALIBS[1][2]; t = (TG..., x_v = 1.0))
gD1 = gs_s(cD1); SD1 = irf_s(cD1, gD1; ρ = ρ_s_cal, H = 12)
P(@sprintf("  D1, s shock, ρ = %.3f, per 1%% rise in s (s̃_0 = 0.01), in %%:", ρ_s_cal))
P(@sprintf("    %4s %9s %9s %9s", "t", "θ̂_t", "û_t", "v̂_t"))
for t in 0:4
    P(@sprintf("    %4d %9.4f %9.4f %9.4f", t, SD1.θh[t+1], SD1.uh[t+1], SD1.vh[t+1]))
end
P("="^100)

text = join(lines, "\n")
println(text)
open(joinpath(@__DIR__, "gs_free_entry_beveridge.txt"), "w") do io
    write(io, text * "\n")
end
