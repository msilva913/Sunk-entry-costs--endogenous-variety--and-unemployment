# R8_test_steady_state.jl
# =============================================================================
# Q1 step 3 test (Notes/Q1_R8_D1_PATHB_plan.md): steady state with partial reposting.
#
# Test 3a. Costless reposting (no Λr_target): the steady state is unchanged from the
#   pre-R8 snapshot (R8_snapshot_pre.jls, case :base).
# Test 3b. Partial reposting, Λ̄_r ∈ {0.95, 0.85, 0.70}: calibrate_shares → steady_state
#   round trip recovers the targets (Λ_r, θ, N, X/Y, δ_e), and every steady-state equation
#   of the draft's app:steady_state holds at 1e-12: the JCC and wage with s𝓡, entry with
#   retired positions, the vacancy LOM with Λ_r, and the resource constraint with X^r.
# Test 3c. The interior closed form for reactivation costs per vacated position,
#   α∫₀^{Q/α} χ dF = ψ_c Q (Λ_r − 1 + p_0), matches numerical quadrature.
#
# Output: R8_test_steady_state.txt
# =============================================================================

include("steady_state.jl")
using Printf, Serialization

tg = merge(TARGETS, (dest_ann = 0.0320, b_ratio = 0.9, x_v = 0.5,
                     p_0 = 0.5, dest_end_frac = 0.5, dest_elast_target = 5.0))

lines = String[]
log!(s) = (push!(lines, s); println(s))
log!("R8 step 3 tests — steady state with partial reposting (D1, PATH B, p_0 = 0.5)")

# ── Test 3a ──────────────────────────────────────────────────────────────────
snap = deserialize("R8_snapshot_pre.jls")[:base].ss
st0  = steady_state(calibrate_shares(tg))
d3a  = maximum(abs(getfield(st0, k) - getfield(snap, k)) for k in (:u, :v, :θ, :N, :δ_e, :x_c, :w, :C))
log!(@sprintf("\nTest 3a (costless): max |SS − pre-R8 SS| over u, v, θ, N, δ_e, x_c, w, C = %.2e  %s",
              d3a, d3a < 1e-10 ? "PASSED" : "FAILED"))

# ── Test 3b / 3c ─────────────────────────────────────────────────────────────
worst = 0.0
for Λr in (0.95, 0.85, 0.70)
    cal = calibrate_shares(merge(tg, (Λr_target = Λr,)))
    st  = steady_state(cal)
    @unpack δ_e, θ, u, v, e, q, f, K, Q, w, w_int, N, Y, Y_c, C, X, X_c, X_r, R_sh, Λ_r = st
    @unpack r, s, ϕ, b, p_0, ψ, f_m, α, κ = cal
    τ   = compute_separation_rate(δ_e, s)
    ψ_c = ψ / (ψ + 1)
    gap = (1 - Λ_r) + ψ_c * (Λ_r - 1 + p_0)
    Qrep = Q * Λ_r - ψ_c * Q * (Λ_r - 1 + p_0)

    res = Dict(
        "Λ_r round trip"  => Λ_r - Λr,
        "θ round trip"    => θ - (tg.f / (1 - δ_e)) / (tg.q / (1 - δ_e)),
        "N normalization" => N - tg.N,
        "X/Y target"      => X / Y - tg.X_Y,
        "δ_e target"      => δ_e - (1 - (1 - tg.dest_ann)^(1 / 12)),
        "JCC (s𝓡)"        => w_int - w - K - s * R_sh - (r + τ) / (1 - δ_e) * (κ + K / q),
        "wage (ϕs𝓡)"      => w - (ϕ * (w_int - K + θ * (K + q * κ) - s * R_sh) + (1 - ϕ) * b),
        "entry eq:e_ss"   => e - (δ_e * (v + 1 - u) + (1 - δ_e) * (1 - Λ_r) * s * (1 - u)),
        "vacancy LOM"     => v - ((1 - δ_e) * ((1 - q) * v + Λ_r * s * (1 - u)) + e),
        "u LOM"           => (1 - δ_e) * f * u - τ * (1 - u),
        "𝓡 definition"    => R_sh - (1 - δ_e) * (Q - Qrep) / (1 + r),
        "RC with X^r"     => Y_c - (C + X + X_c + X_r),
        "Λ_r = F(Q/α)"    => Λ_r - (1 - p_0 + p_0 * (Q / (α * f_m))^ψ),
        "w normalization" => w - tg.w,
    )
    # 3c: reactivation cost per vacated position by quadrature of α∫₀^{Q/α} χ dF
    #     (dF = p_0 ψ χ^(ψ-1)/χ_m^ψ dχ on the power-law part; the atom at 0 costs nothing)
    xbar = Q / α
    n    = 200_000
    grid = range(0, xbar; length = n + 1)
    dens(χ) = χ == 0 ? 0.0 : p_0 * ψ * χ^(ψ - 1) / f_m^ψ
    # χ·dens(χ) = p_0 ψ χ^ψ / χ_m^ψ is continuous at 0, so the trapezoid rule is accurate
    integrand = [χ * dens(χ) for χ in grid]
    quad = α * (xbar / n) * (sum(integrand) - (integrand[1] + integrand[end]) / 2)
    res["cost integral (quad vs closed)"] = (quad - ψ_c * Q * (Λ_r - 1 + p_0)) / (ψ_c * Q)

    log!(@sprintf("\nΛ̄_r = %.2f:  α = %.5g   Q = %.5f   𝓡/Q = %.4f   X^r/Y = %.4f%%   ϕ = %.4f   b/w_int = %.4f",
                  Λr, α, Q, R_sh / Q, 100 * X_r / Y, ϕ, b / w_int))
    for k in sort(collect(keys(res)))
        tol = startswith(k, "cost integral") ? 1e-6 : 1e-10
        ok  = abs(res[k]) < tol
        k == "cost integral (quad vs closed)" || (global worst = max(worst, abs(res[k])))
        log!(@sprintf("    %-32s %10.2e  %s", k, res[k], ok ? "ok" : "FAIL"))
    end
end
log!(@sprintf("\nTest 3b %s (worst identity residual %.2e, tolerance 1e-10)", worst < 1e-10 ? "PASSED" : "FAILED", worst))
write("R8_test_steady_state.txt", join(lines, "\n") * "\n")
