# double_limit_check.jl
# =============================================================================
# Check of prop:double_limit claim (i) (app:limiting_cases). DS-CES, p_0 = 0, f_e = f̄_e/ε.
#
# Exact steady state of eq:N_lom_eq and eq:N_euler_eq (Λ = 1, X^c = 0, β = 1/(1+r)):
#   Euler:  Y^c/N = f_e ρ(N)(r+δ)/[(1−δ)(μ−1)]
#   LOM:    δN/(1−δ) = N^e = zL/f_e − Y^c/(f_e ρ(N))
#   ⇒  N = (μ−1) z L (1−δ) / [f_e (r + δμ)]          (= eq:curve_res_cons; no ρ(N))
#   ε → ∞ with f_e = f̄_e/ε:  N → N_lim ≡ z L (1−δ) / [f̄_e (r+δ)]
#
# The draft's eq:N_ss_limit, δN = zL lnN/f̄_e, comes from replacing Y^c by its limit zL
# inside N^e = (ε/f̄_e)(zL − Y^c/ρ): an O(1/ε) error times ε. Its roots are printed for
# comparison; it has 0 or 2 roots, never exactly one.
#
# OUTPUT  double_limit_check.txt
# =============================================================================

using Printf

const r = (1.04)^(1 / 12) - 1        # monthly, as in beveridge_free_entry.jl
const δ = 1 - (1 - 0.0320)^(1 / 12)  # D1
const L = 0.93;  const z = 1.0
N_exact(ε, f̄) = (1 / (ε - 1)) * z * L * (1 - δ) / ((f̄ / ε) * (r + δ * ε / (ε - 1)))
N_lim(f̄) = z * L * (1 - δ) / (f̄ * (r + δ))

"Exact steady-state residuals of eq:N_lom_eq and eq:N_euler_eq at N (σ = 0)."
function resid(N, ε, f̄)
    fe = f̄ / ε;  μ = ε / (ε - 1);  ρN = N^(1 / (ε - 1));  β = 1 / (1 + r)
    Yc = fe * ρN * N * (r + δ) / ((1 - δ) * (μ - 1))           # from the Euler equation
    lom = N - (1 - δ) * (N + z * L / fe - Yc / (fe * ρN))
    eul = fe * ρN - β * (1 - δ) * (fe * ρN + Yc / N * (μ - 1))
    return lom, eul
end

"Roots of the draft's δN = c·lnN, c = zL/f̄ (0 or 2 roots)."
function draft_roots(f̄)
    c = z * L / f̄;  g(N) = c * log(N) - δ * N
    Nmax = c / δ                                                # argmax of g
    g(Nmax) < 0 && return Float64[]
    bis(lo, hi) = (for _ in 1:200; m = (lo + hi) / 2; sign(g(m)) == sign(g(lo)) ? (lo = m) : (hi = m); end; (lo + hi) / 2)
    return [bis(1.0, Nmax), bis(Nmax, 1e12)]
end

lines = String[];  P(s) = push!(lines, s)
P("="^90)
P(@sprintf("prop:double_limit (i): steady-state N under f_e = f̄_e/ε  (r = %.5f, δ = %.5f, L = %.2f, z = 1)", r, δ, L))
P("="^90)
for f̄ in (z * L * (1 - δ) / (r + δ), 0.05, 0.5)                # first: f̄ normalizing N_lim = 1
    P("")
    P(@sprintf("f̄_e = %.4f:  N_lim = zL(1−δ)/(f̄_e(r+δ)) = %.4f", f̄, N_lim(f̄)))
    P(@sprintf("  %10s %14s %14s %12s %12s", "ε", "N exact", "N/N_lim − 1", "LOM resid", "Euler resid"))
    for ε in (5.0, 10.0, 100.0, 1e3, 1e5)
        N = N_exact(ε, f̄);  lr, er = resid(N, ε, f̄)
        @assert abs(lr) < 1e-10 * z * L * ε / f̄ && abs(er) < 1e-10 * f̄ / ε "exact SS fails at ε = $ε"   # LOM cancels terms of order zLε/f̄
        P(@sprintf("  %10g %14.6f %14.2e %12.1e %12.1e", ε, N, N / N_lim(f̄) - 1, lr, er))
    end
    @assert isapprox(N_exact(1e7, f̄), N_lim(f̄); rtol = 1e-5)
    rts = draft_roots(f̄)
    P("  draft eq:N_ss_limit δN = zL lnN/f̄_e: " *
      (isempty(rts) ? "no root" : @sprintf("two roots, N = %.4f and %.4f", rts...)))
    for N in rts                                                 # draft roots are not steady states
        lr, _ = resid(N, 1e5, f̄)
        P(@sprintf("    at N = %.4f the exact LOM residual (ε = 1e5) is %.3e  (vs N scale %.3e)", N, lr, N))
    end
end
P("="^90)

text = join(lines, "\n");  println(text)
open(joinpath(@__DIR__, "double_limit_check.txt"), "w") do io; write(io, text * "\n"); end
