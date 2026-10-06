# gs_free_entry_beveridge.jl
# =============================================================================
# δ shock at free entry (ξ → ∞) under the conditions of prop:independence
# (σ = 0, ζ = 0, p_0 = 0): the labor-market block closes on its own and is the
# Gabrovski-Silva economy with a time-varying δ_t. Draft timing (S10): δ_t is known
# at the start of t and acts at the end of t, so u_t and v_pre,t are predetermined
# and u_{t+1}, v_pre,t+1 carry the destruction. Companion to
# Notes/beveridge_free_entry_GS_proposition.md.
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
    return (; a, Φ, λ, c_u, ω, κu, v1, slope, uc, Γ1, ΓP)
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
P("="^100)

text = join(lines, "\n")
println(text)
open(joinpath(@__DIR__, "gs_free_entry_beveridge.txt"), "w") do io
    write(io, text * "\n")
end
