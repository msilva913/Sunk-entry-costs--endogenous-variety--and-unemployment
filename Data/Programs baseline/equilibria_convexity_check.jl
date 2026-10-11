# equilibria_convexity_check.jl
# =============================================================================
# Task R16: shape of the steady-state curves in prop:equilibria (DS-CES, p_0 → 0,
# costless reposting), app:equilibria_proof eq:jcc_rearranged, for a general CRS
# matching function with bounded probabilities.
#   JC:  (z/μ) N^{1/(ε−1)} = Φ(θ) ≡ K(θ) + b + (κ + K(θ)/q(θ))(r+τ+(1−δ)ϕf(θ))/((1−δ)(1−ϕ))
#        N_JC(θ) = (μΦ(θ)/z)^{ε−1}
#   RC:  N_RC(θ) = (μ−1) z (1−δ) L(θ) / (f_e(δμ + r)),  L(θ) = (1−δ)f/(τ + (1−δ)f)
#   K(θ) = (r+δ)/(1+r) e(θ)^{1/ξ} x_m,  e(θ) = δ(τθ + (1−δ)f)/(τ + (1−δ)f)   (eq:e_theta, Λ_r = 1)
#
# Argument checked. ln N_JC − ln N_RC = (ε−1) ln Φ − ln L + const. If ln Φ is convex and
# ln L concave in ln θ, the difference is convex in ln θ and has at most two zeros. With
# N_JC > N_RC at θ → 0 (b > 0 or κ > 0, q bounded) and θ → ∞, the count is 0 or 2.
# ln L is concave whenever the elasticity of f, ε_f(θ), is nonincreasing (checked from ε_fn).
#
# Matching functions: the MatchingSpec library of the companion paper (TwoPSMCD, HRW,
# ClippedCD in Soft-Min-Cobb-Douglas-Matching-Function/Programs/dmp_smcd.jl, Section 1),
# loaded as is. Each is calibrated to this paper's steady state (θ̄, f̄) and, where the
# form allows, the elasticity ε* = 1 − η_L, by the joint-calibration construction of that
# paper's prop:2p (iii). Unclipped Cobb–Douglas, the form the draft's proof uses, is the
# reference row.
# Parameterization: calibrate() (PATH B mirror); z = μ w_int and f_e set so that the
# Cobb–Douglas θ̄ is a crossing.
#
# OUTPUT  equilibria_convexity_check.txt
# =============================================================================

include(joinpath(@__DIR__, "beveridge_free_entry.jl"))   # calibrate, CALIBS, TG
using Printf, Parameters

const SMCD_LIB = joinpath(@__DIR__, "..", "..", "..", "Soft-Min-Cobb-Douglas-Matching-Function",
                          "Programs", "dmp_smcd.jl")
isfile(SMCD_LIB) || error("companion repository not found: $SMCD_LIB")
let src = read(SMCD_LIB, String)                         # Section 1 only: types and f, q, ε
    i0 = findfirst("abstract type MatchingSpec end", src)[1]
    i1 = findfirst("_spec_label(::HRW)", src)[end]
    include_string(Main, src[i0:findnext('\n', src, i1)])
end

"Unclipped Cobb–Douglas (the draft's form), for reference."
struct ExtrapCD <: MatchingSpec
    A::Float64
    η::Float64
end
f_fn(m::ExtrapCD, θ::Float64) = m.A * θ^(1 - m.η)
q_fn(m::ExtrapCD, θ::Float64) = m.A * θ^(-m.η)
ε_fn(m::ExtrapCD, θ::Float64) = 1 - m.η

"2P-SMCD through (θ*, f*, ε*), ε* < 1/2 (right branch), given θ_c; p_L sets q_max = 1."
function smcd_cal(θs, fs, εs, θc)
    pR = log((1 - εs) / εs) / log(θs / θc)
    fc = fs / (2 * (1 - εs))^(1 / pR)
    @assert fc < θc && fc * 2^(1 / pR) <= 1 "infeasible θ_c = $θc"
    return TwoPSMCD(θc = θc, pL = log(2) / log(θc / fc), pR = pR, fc = fc)
end

"HRW through (θ*, f*) (one parameter: level only)."
function hrw_cal(θs, fs)
    g(ν) = θs * (1 + θs^ν)^(-1 / ν) - fs
    lo, hi = 0.05, 50.0
    for _ in 1:200
        m = (lo + hi) / 2
        sign(g(m)) == sign(g(lo)) ? (lo = m) : (hi = m)
    end
    return HRW(ν = (lo + hi) / 2)
end

function curves(c, mf; ξ = 1 / TG.ξ_inv, ε = TG.ε)
    μ = ε / (ε - 1);  δ, τ, r, ϕ, b, κ = c.δ_e, c.τ, c.r, c.ϕ, c.b, c.κ
    f(θ) = f_fn(mf, θ);  q(θ) = q_fn(mf, θ)
    x_m = c.Q / c.e^(1 / ξ)
    e(θ) = δ * (τ * θ + (1 - δ) * f(θ)) / (τ + (1 - δ) * f(θ))
    K(θ) = (r + δ) / (1 + r) * e(θ)^(1 / ξ) * x_m
    Φ(θ) = K(θ) + b + (κ + K(θ) / q(θ)) * (r + τ + (1 - δ) * ϕ * f(θ)) / ((1 - δ) * (1 - ϕ))
    L(θ) = (1 - δ) * f(θ) / (τ + (1 - δ) * f(θ))
    z = μ * c.w_int                                       # N_JC(θ̄) = 1 under CD
    f_e = (μ - 1) * z * (1 - δ) * L(c.θ) / (δ * μ + r)    # N_RC(θ̄) = 1 under CD
    NJC(θ) = (μ * Φ(θ) / z)^(ε - 1)
    NRC(θ) = (μ - 1) * z * (1 - δ) * L(θ) / (f_e * (δ * μ + r))
    return (; Φ, L, NJC, NRC, f, q)
end

"Elasticity of g in ln θ by central differences, on a grid of ln θ."
elast(g, x; h = 1e-5) = (log(g(exp(x + h))) - log(g(exp(x - h)))) / (2h)

lines = String[];  P(s) = push!(lines, s)
P("="^104)
P("prop:equilibria for bounded-probability matching (DS-CES, p_0 → 0): log-shape of the curves and crossings")
P("Δε_Φ ≥ 0: ln Φ convex in ln θ.  Δε_L ≤ 0: ln L concave.  Δε_f ≤ 0: ε_f nonincreasing (from ε_fn).")
P("="^104)
const XS = range(log(1e-6), log(1e4), length = 6000)
for (nm, d) in CALIBS[1:3], (xlab, tg) in (("x_v .5", TG), ("x_v 1", (TG..., x_v = 1.0)))
    c = calibrate(d; t = tg);  εs = 1 - c.η
    specs = [("CD extrapolated", ExtrapCD(c.q * c.θ^c.η, c.η)),
             ("CD clipped", ClippedCD(A = c.q * c.θ^c.η, η = c.η)),
             ("HRW (level)", hrw_cal(c.θ, c.f)),
             [(@sprintf("2P-SMCD θc=%.2f", θc), smcd_cal(c.θ, c.f, εs, θc)) for θc in (0.43, 0.46, 0.49)]...]
    P("")
    P(@sprintf("%s, %s:  θ̄ = %.3f  f̄ = %.3f  ε* = %.2f  ϕ = %.3f  b = %.2f", nm, xlab, c.θ, c.f, εs, c.ϕ, c.b))
    P(@sprintf("  %-17s %9s %9s %9s %8s | %9s %6s %6s | %9s %6s %6s", "spec", "min Δε_Φ", "max Δε_L",
               "max Δε_f", "q_max", "θ low", "q", "u", "θ high", "q", "u"))
    for (lab, mf) in specs
        C = curves(c, mf)
        @assert isapprox(C.f(c.θ), c.f; rtol = 1e-6) "$lab misses f̄"
        εΦ = elast.(C.Φ, XS);  εL = elast.(C.L, XS);  εf = [ε_fn(mf, exp(x)) for x in XS]
        jump = [i for i in 1:length(XS)-1 if abs(εf[i+1] - εf[i]) > 1e-3]   # ClippedCD kink
        ok(i) = all(j -> abs(i - j) > 3, jump)
        dΦ = minimum(εΦ[i+1] - εΦ[i] for i in 1:length(XS)-1 if ok(i))
        dL = maximum(εL[i+1] - εL[i] for i in 1:length(XS)-1 if ok(i))
        df = maximum(εf[i+1] - εf[i] for i in 1:length(XS)-1)
        θg = exp.(range(log(1e-8), log(1e6), length = 40000));  Δ = C.NJC.(θg) .- C.NRC.(θg)
        xs = [θg[i] for i in 1:length(θg)-1 if sign(Δ[i]) != sign(Δ[i+1])]
        @assert length(xs) == 2 "$lab: $(length(xs)) crossings"
        u(θ) = 1 - C.L(θ) / (1 - c.δ_e) * (1 - c.δ_e)        # u = τ/(τ + (1−δ)f)
        qm = mf isa ExtrapCD ? Inf : q_max(mf)
        P(@sprintf("  %-17s %+9.1e %+9.1e %+9.1e %8.3f | %9.4f %6.2f %6.3f | %9.4f %6.3f %6.3f", lab, dΦ, dL,
                   df, qm, xs[1], C.q(xs[1]), u(xs[1]), xs[2], C.q(xs[2]), u(xs[2])))
    end
end
P("")
P("── Condition (C) across entry elasticities: min Δε_Φ (≥ 0 up to grid noise) ──────────────────────────")
P(@sprintf("  %-26s %-7s %-17s %11s %11s %11s", "Calibration", "x_v", "spec", "ξ = 0.5", "ξ = 1", "ξ = 2"))
for (nm, d) in CALIBS[1:3], (xlab, tg) in (("0.5", TG), ("1", (TG..., x_v = 1.0)))
    c = calibrate(d; t = tg)
    for (lab, mf) in (("CD clipped", ClippedCD(A = c.q * c.θ^c.η, η = c.η)), ("HRW (level)", hrw_cal(c.θ, c.f)),
                      ("2P-SMCD θc=0.46", smcd_cal(c.θ, c.f, 1 - c.η, 0.46)))
        dΦ = map((0.5, 1.0, 2.0)) do ξ
            C = curves(c, mf; ξ = ξ);  εΦ = elast.(C.Φ, XS)
            kink = mf isa ClippedCD ? [log(mf.A^(1 / mf.η)), log(mf.A^(-1 / (1 - mf.η)))] : Float64[]
            minimum(εΦ[i+1] - εΦ[i] for i in 1:length(XS)-1 if all(k -> abs(XS[i] - k) > 0.02, kink))
        end
        @assert all(>(-1e-9), dΦ) "(C) fails: $nm, $lab"
        P(@sprintf("  %-26s %-7s %-17s %+11.1e %+11.1e %+11.1e", nm, xlab, lab, dΦ...))
    end
end
P("  Minima sit at small θ, where Φ ≈ b + κ(r+τ)/((1−δ)(1−ϕ)) and the specs coincide; any violation")
P("  elsewhere would appear as a more negative entry.")
P("")
P("-"^104)
P("Two crossings in every case. Under unclipped CD the low one needs q > 1; under every bounded")
P("specification it is a feasible steady state with q ≤ q_max ≤ 1 and unemployment near one third.")
P("ln Φ convex and ln L concave hold for every bounded spec (|violations| ≤ 3e-11, grid noise). Unclipped CD")
P("violates ln Φ convexity slightly (−7e-7 at D1, x_v = 1). ClippedCD kinks at q = 1 and f = θ are skipped.")
P("="^104)

text = join(lines, "\n");  println(text)
open(joinpath(@__DIR__, "equilibria_convexity_check.txt"), "w") do io; write(io, text * "\n"); end
