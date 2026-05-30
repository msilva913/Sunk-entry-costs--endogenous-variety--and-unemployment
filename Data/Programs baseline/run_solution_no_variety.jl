# run_solution_no_variety.jl
# =============================================================================
# Mechanism comparison C: role of VARIETY (love-of-variety / price-index effect)
#
# Compares the full baseline against a counterfactual in which the variety
# externality is suppressed: ρ_t ≡ 1 (no taste for diversity, ζ = 0).
#
#   (1) Full baseline
#         ρ_t = N_t^{1/(ε−1)} — standard DS-CES relative price
#         Same as baseline in Comparisons A and B
#
#   (2) No-variety counterfactual (ζ = 0)
#         ρ_t ≡ 1  — relative price fixed, no love-of-variety externality
#         Markup μ = ε/(ε−1) is UNCHANGED (same ε = 4.3)
#         Consequence: ν_f = f_e/μ (constant), w_int = z/μ (no N amplification)
#         Free-entry condition ν_f = f_e/μ still pins f_e, but f_e adjusts to
#         compensate for the removed ρ(N) factor → f_e shifts, N̄ need not equal 1
#
# COMPARISON DESIGN:
#   Identification principle: hold fixed markup, τ (total separation rate), and
#   all other structural parameters. Shut off only the N → ρ → w_int channel.
#
#   (a) Same τ = 0.031: keeps steady-state u, v, θ approximately equal, so
#       pp-deviation IRF comparisons are valid.
#   (b) Same ε = 4.3, same b_ratio, x_v: only structural change is ρ ≡ 1.
#   (c) f_e adjusts (via Stage 4 of calibration) to be consistent with ρ = 1
#       and the free-entry condition. N̄ may differ from 1 — not targeted.
#
# Implementation:
#   - After gen_model_equations(), replace f[13] (ρ = N^{1/(ε-1)}) with ρ = 1.
#   - In calibration, set ρ_val = 1.0 instead of N^{1/(ε-1)} throughout Stage 4.
#   - All other equations (JCC, Nash wage, Euler) use ρ as-is; with ρ=1 they
#     mechanically suppress the N-feedback on w_int, ν_f, and χ^c.
#
# IRFs computed for z (technology) and δ (destruction) shocks.
# Serialized output for plotting in plot_no_variety_comparison.jl.
# =============================================================================

include("run_solution_core.jl")

# ── Shared shock process parameters (monthly, from part6b) ────────────────────
ρ_z = 0.902
σ_z = 0.0092
ρ_δ = 0.592
σ_δ = 0.0669
ρ_s = 0.90
σ_s = 0.010

T_IR        = 60    # 5 years at monthly frequency
flag_IR     = true
flag_logdev = true

# =============================================================================
# Helper: calibrate_shares_no_variety
# Wraps calibrate_shares logic with ρ ≡ 1 substituted throughout Stage 4.
# All other stages (determining δ_e, s, π_s, ψ, cons) are identical.
# Only Stage 4 differs: ρ = 1 replaces N^{1/(ε-1)}, so z = μ * w_int and
# ν_f = f_e/μ. The free-entry condition still pins f_e analytically.
# =============================================================================
function calibrate_shares_no_variety(targets)
    @unpack X_Y, dest_ann, dest_end_frac, p_0, f, η_L, q, sep, b_ratio,
            x_v, ξ_inv, ε, r_ann, σ, N, w = targets
    Xc_Y              = get(targets, :Xc_Y, 0.10)
    dest_elast_target = get(targets, :dest_elast_target, nothing)
    use_dest_elast    = !isnothing(dest_elast_target) && p_0 > 0.0

    # Stage 1 (identical to baseline)
    δ_e      = 1 - (1 - dest_ann)^(1 / 12)
    δ        = (1 - dest_end_frac) * δ_e
    surv_prob = (1 - δ_e) / (1 - δ)

    μ  = compute_markup(ε)
    τ  = sep
    s  = (τ - δ_e) / (1 - δ_e)
    r  = (1 + r_ann)^(1 / 12) - 1

    f_corr = f / (1 - δ_e)
    q_corr = q / (1 - δ_e)
    θ  = f_corr / q_corr
    u  = τ / (τ + (1 - δ_e) * f_corr)
    v  = θ * u
    L  = 1 - u
    A  = f_corr / θ^(1 - η_L)

    e   = δ_e * (v + 1 - u)
    N_e_share = δ_e / (1 - δ_e)   # N_e / N ratio; N itself is an outcome here
    b   = b_ratio * w

    # ── KEY CHANGE: ρ ≡ 1 (no variety externality) ────────────────────────────
    ρ = 1.0

    ζ_ss = (p_0 > 0.0) ? (surv_prob - (1 - p_0)) / p_0 : 0.0

    # Stages 2-3 (identical to baseline — π_s, cons, ψ)
    if p_0 == 0.0
        cons = 0.0
        ψ    = 1.5
        π_s  = 1 / ε
    elseif use_dest_elast
        0 < ζ_ss < 1 || error("ζ_ss = $(round(ζ_ss,digits=6)) not in (0,1).")
        ψ    = dest_elast_target * (1 - ζ_ss) / ζ_ss
        cons = ψ / (1 + ψ) * p_0
        π_s  = (μ - 1) / μ * (1 - cons) * (r + δ_e) / (r + δ_e + cons * (1 - δ_e))
        Xc_Yc_implied = 1 / ε - π_s
        π_s ≤ 0.0 && error("dest_elast_target infeasible: π_s = $(round(π_s,digits=4)) ≤ 0.")
        println("  [NO-VAR PATH B] ψ=$(round(ψ,digits=4)), π_s=$(round(π_s,digits=4))")
    else
        if Xc_Y == 0.0
            Xc_Yc = 0.0
            π_s   = 1 / ε
        else
            function loss_xc(xc_yc)
                π_s_trial = 1 / ε - xc_yc
                Yc_YG = (r + δ_e) / (r + δ_e + δ_e * π_s_trial)
                YG_Y  = 1 + X_Y + Xc_Y
                return Xc_Y / (Yc_YG * YG_Y) - xc_yc
            end
            Xc_Yc = find_zero(loss_xc, (0.01, 0.5))
            π_s   = 1 / ε - Xc_Yc
        end
        function loss_psi(cons_trial)
            return (μ - 1) / μ * (1 - cons_trial) * (r + δ_e) /
                   (r + δ_e + cons_trial * (1 - δ_e)) - π_s
        end
        cons = find_zero(loss_psi, 0.1)
        ψ_c  = cons / p_0
        ψ    = ψ_c / (1 - ψ_c)
    end

    L_c = (r + δ_e) * L / (r + δ_e + δ_e * π_s * μ)
    L_e = L - L_c
    surplus_ratio = (r + τ) / (1 - δ_e) / (q_corr * x_v)

    # Stage 4: Solve for Q (ρ ≡ 1 throughout)
    # With ρ = 1: z = μ * w_int; ν_f = f_e/μ (no N dependence)
    # N is now a derived object (from free entry), not pinned to 1
    function loss_Q(Q_trial)
        Q_trial = abs(Q_trial)
        K     = Q_trial * (r + δ_e) / (1 + r)
        κ     = (1 - x_v) / x_v * K / q_corr
        X     = e / (1 + ξ_inv) * Q_trial + κ * q_corr * v
        w_int = surplus_ratio * K + w + K
        z     = (μ / ρ) * w_int                    # ρ = 1 → z = μ * w_int
        # N from free-entry: ν_f = ρ*f_e/μ = f_e/μ, and from profit:
        #   f_e = π_s * z * L_c * (1-δ_e) * μ / (N * (r+δ_e))
        # With ρ=1, the calibrated f_e depends on an assumed N. We set N_cal
        # from the resource constraint using N = 1 as a normalization for f_e,
        # then recover implied N via free entry ex post (it is an outcome).
        # Alternatively, use the same formula as baseline with N=1 for f_e:
        f_e   = π_s * z * L_c * (1 - δ_e) * μ / (1.0 * (r + δ_e))
        ν_f   = ρ * f_e / μ                        # = f_e/μ (constant)
        d_f   = (r + δ_e) / (1 - δ_e) * ν_f
        Y_c   = ρ * z * L_c                        # = z * L_c
        # N_e = N * δ_e/(1-δ_e); N implicit. Use N=1 for x_c and X_c computation:
        N_cal = 1.0
        x_c   = Y_c / (ε * N_cal) + ν_f
        X_c   = N_cal * cons * x_c
        C     = Y_c - X - X_c
        N_e_cal = δ_e / (1 - δ_e) * N_cal
        Y_new = C + ν_f * N_e_cal
        Y_rec = X / X_Y
        loss_val = 100 * (Y_rec - Y_new) / (Y_rec + Y_new)
        return loss_val, (; w_int, κ, z, f_e, K, d_f, ν_f, x_c, X_c, X, C, Y_c, Y=Y_new)
    end

    Q_opt = find_zero(Q -> loss_Q(Q)[1], 1.0)
    Q     = abs(Q_opt)
    _, Q_out = loss_Q(Q)
    @unpack w_int, κ, z, f_e, K, d_f, ν_f, x_c, X_c, X, C, Y_c = Q_out

    # Stage 5: Recover distribution parameters (same as baseline)
    if p_0 == 0.0
        ζ       = 0.0
        dest_el = 0.0
        f_m     = 1.0
    else
        ζ       = (surv_prob - (1 - p_0)) / p_0
        dest_el = ψ * ζ / (1 - ζ)
        f_m     = x_c / ζ^(1 / ψ)
    end

    ϕ   = (w - b) / (w_int - K + θ * (K + q_corr * κ) - b)
    x_m = Q / e^ξ_inv

    return (;
        f_e, δ, z, b, ϕ, r, σ, ε, A, η_L,
        κ, ξ_inv, x_m, s, ψ, p_0, f_m, dest_el
    )
end

# =============================================================================
# Helper: build model, solve, compute IRFs for z and δ shocks
# For the no-variety arm, f[13] is replaced with ρ = 1 after equation generation.
# =============================================================================
function solve_and_irf_C(targets_variant; no_variety=false)
    if no_variety
        cal = calibrate_shares_no_variety(targets_variant)
    else
        cal = calibrate_shares(targets_variant)
    end
    @unpack f_e, δ, s, z, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, x_m, ψ, f_m, p_0 = cal

    if s < -1e-8
        error("Negative s = $s: τ < δ_e, infeasible calibration.")
    end
    s = max(s, 0.0)

    zbar = z; δbar = δ; sbar = s

    PAR = [f_e; zbar; δbar; sbar; b; ϕ; r; σ; ε; A; η_L; κ; ξ_inv; x_m; ψ; f_m; p_0;
           ρ_z; σ_z; ρ_δ; σ_δ; ρ_s; σ_s]

    SS_symbolic = SS_symbolics(parameters, targets_variant)
    SS_numeric_val = SS_numeric(PAR, targets_variant)

    # Generate model equations; override f[13] for no-variety arm
    f_eqs = gen_model_equations()
    if no_variety
        # Replace f[13]: ρ - N^{1/(ε-1)} = 0  →  ρ - 1 = 0
        # ρ is in the control vector y; ε is a symbolic parameter.
        # After substitution ρ = 1 for all N, ρ is still a control variable
        # but pinned to 1 by this equation rather than CES aggregation.
        f_eqs[13] = ρ - 1
    end

    model = (parameters = parameters, estimate = estimate, estimation = position,
             npar      = length(parameters), ns = length(estimate),
             priors    = priors,
             x = x, y = y, xp = xp, yp = yp, variables = variables,
             varnames  = varnames,
             nx = nx, ny = ny, nvar = nvar,
             e = ex, eta = eta,
             ne = ne,
             f = f_eqs,
             nf = nvar,
             SS = SS_symbolic, PAR_SS = parameters[:],
             flag_order = flag_order, flag_deviation = flag_deviation,
             flag_SSsolver = flag_SSsolver)
    process_model(model)

    sol = solution_interface(model, PAR, SS_numeric_val)
    @unpack ss, SS, sol_mat = sol

    eta_z_col = reshape([0.0, 0.0, 0.0, σ_z, 0.0, 0.0], nx, 1)
    eta_δ_col = reshape([0.0, 0.0, 0.0, 0.0, σ_δ, 0.0], nx, 1)

    model_1shock = (; model..., ne = 1)

    irf_z = simulate_model(model_1shock, sol_mat, T_IR, eta_z_col, SS, flag_IR, flag_logdev)
    irf_z = 100 .* DataFrame(irf_z, varnames)

    irf_δ = simulate_model(model_1shock, sol_mat, T_IR, eta_δ_col, SS, flag_IR, flag_logdev)
    irf_δ = 100 .* DataFrame(irf_δ, varnames)

    # Exit flow (δ_e × N in log dev) — kept for consistency with B
    irf_z[!, :exit_flow] = irf_z.δ_e .+ irf_z.N
    irf_δ[!, :exit_flow] = irf_δ.δ_e .+ irf_δ.N

    return (ss=ss, irf_z=irf_z, irf_δ=irf_δ, cal=cal)
end

# =============================================================================
# (1) Full baseline — identical targets to Comparisons A and B
# (same as baseline in run_solution_endog_exit.jl: endog exit active,
#  dest_elast_target=5, b_ratio=0.9, x_v=0.5)
# =============================================================================
println("\n" * "="^60)
println("Solving BASELINE model (full variety, ρ = N^{1/(ε-1)})")
println("="^60)

targets_base = (TARGETS...,
    dest_end_frac     = 0.5,
    p_0               = 0.5,
    dest_elast_target = 5.0,
    b_ratio           = 0.9,
    x_v               = 0.5,
)

out_base = solve_and_irf_C(targets_base; no_variety=false)
println("Baseline SS: u=$(round(out_base.ss.u, digits=4)), " *
        "δ_e=$(round(out_base.ss.δ_e, digits=5)), " *
        "ρ=$(round(out_base.ss.ρ, digits=4)), " *
        "N=$(round(out_base.ss.N, digits=4))")

# =============================================================================
# (2) No-variety counterfactual — ρ ≡ 1, same τ, same ε (markup preserved)
# =============================================================================
println("\n" * "="^60)
println("Solving NO-VARIETY model (ρ ≡ 1, ζ = 0)")
println("="^60)

# Same targets as baseline, only structural change is ρ ≡ 1 in equations.
# dest_elast_target=5 is kept so the endogenous exit margin is active with
# the same elasticity — isolating only the variety channel.
targets_novar = (TARGETS...,
    dest_end_frac     = 0.5,
    p_0               = 0.5,
    dest_elast_target = 5.0,
    b_ratio           = 0.9,
    x_v               = 0.5,
)

out_novar = solve_and_irf_C(targets_novar; no_variety=true)
println("No-variety SS: u=$(round(out_novar.ss.u, digits=4)), " *
        "δ_e=$(round(out_novar.ss.δ_e, digits=5)), " *
        "ρ=$(round(out_novar.ss.ρ, digits=4)), " *
        "N=$(round(out_novar.ss.N, digits=4))")

# =============================================================================
# Serialize for plotting
# =============================================================================
output_C = (
    base  = (ss=out_base.ss,  irf_z=out_base.irf_z,  irf_δ=out_base.irf_δ,
             cal=out_base.cal,  targets=targets_base),
    novar = (ss=out_novar.ss, irf_z=out_novar.irf_z, irf_δ=out_novar.irf_δ,
             cal=out_novar.cal, targets=targets_novar),
)
serialize("irf_no_variety.jls", output_C)
println("\nSaved: irf_no_variety.jls")

# =============================================================================
# Quick diagnostic
# =============================================================================
println("\n── Steady-State Comparison ──────────────────────────────────")
for field in [:u, :v, :θ, :δ_e, :N, :N_e, :ν_f, :x_c, :ρ, :w_int]
    b_val  = getfield(out_base.ss,  field)
    nv_val = getfield(out_novar.ss, field)
    println("  $field:  baseline=$(round(b_val, digits=5))  " *
            "no_var=$(round(nv_val, digits=5))")
end
println("  s:  baseline=$(round(out_base.cal.s, digits=5))  " *
        "no_var=$(round(out_novar.cal.s, digits=5))")

# =============================================================================
# NOTES ON COMPARISON C DESIGN
# =============================================================================
#
# ── What ρ ≡ 1 shuts off ─────────────────────────────────────────────────────
#
# Standard DS-CES: ρ_t = N_t^{1/(ε−1)} > 1 when N_t > 1.
# This means:
#   (1) N → ρ → w_int: larger variety raises the real marginal revenue product
#       of labor → higher wages → tighter JCC → amplified θ response to shocks.
#   (2) N → ρ → ν_f = ρ·f_e/μ: firm value rises with N, making entry forward-
#       looking and entry responses partially self-reinforcing.
#   (3) N → ρ → χ^c: the threshold rises with ρ, dampening endogenous exit when N
#       is high (stabilizing feedback).
#
# With ρ ≡ 1: w_int = z/μ (no N dependence), ν_f = f_e/μ (constant).
# Channels (1)–(3) are all shut off simultaneously. The IRF comparison
# therefore shows the net variety amplification: difference in u and v responses
# is attributable to the aggregate N → ρ externality.
#
# ── Comparison to A and B ────────────────────────────────────────────────────
#
# Comparison A (high δ_e):  isolates the role of the destruction rate level
#   via the vacancy LOM survival coefficient (1−δ_e) and asset duration D(δ_e,ρ).
# Comparison B (exog exit):  isolates the x_c amplification channel.
# Comparison C (ρ ≡ 1):      isolates the N → ρ variety externality channel.
#
# Together A+B+C decompose the three main amplification mechanisms of the model.
# =============================================================================
