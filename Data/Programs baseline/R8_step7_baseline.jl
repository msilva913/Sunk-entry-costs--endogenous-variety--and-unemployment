# R8_step7_baseline.jl
# =============================================================================
# Q1 step 7 (Notes/Q1_R8_D1_PATHB_plan.md). Exploratory, NOT draft-bound (§5.3 is
# regenerated in Q5 after the Q2 VAR).
#
# Part A. Baseline at the code defaults (D1, PATH B, part6b shocks) over the reactivation
#   grid Λ̄_r ∈ {costless, 0.95, 0.85, 0.70}: impact and month-12 responses of u and v to
#   1-SD z, δ, s shocks; simulated cor(u, v) with all shocks and shock by shock
#   (monthly → quarterly averages → HP λ = 1,600 on log deviations, as simulated_moments).
#   Data: cor(u, v) = −0.804 (context/findings.md).
# Test 6 (the prop:fe limit) is in R8_test6_free_entry.jl.
#
# Output: R8_step7_baseline.txt
# =============================================================================

include("run_solution_core.jl")
using Printf, Statistics

const T_IR = 60
const T_SM = 60_000
const LAM  = 1_600.0
const GRID = (nothing, 0.95, 0.85, 0.70)

lines = String[]
log!(s) = (push!(lines, s); println(s))

function build_model(tg; no_var = false)
    f = gen_model_equations()
    no_var && (f[13] = ρ - 1)
    model = (parameters = parameters, estimate = estimate, estimation = position,
             npar = length(parameters), ns = length(estimate), priors = priors,
             x = x, y = y, xp = xp, yp = yp, variables = variables, varnames = varnames,
             nx = nx, ny = ny, nvar = nvar, e = ex, eta = eta, ne = ne,
             f = f, nf = nvar, SS = SS_symbolics(parameters, tg), PAR_SS = parameters[:],
             flag_order = flag_order, flag_deviation = flag_deviation,
             flag_SSsolver = flag_SSsolver)
    process_model(model)
    return model
end

function solve(model, tg, shk = SHOCKS_BASE)
    cal = calibrate_shares(tg)
    PAR = build_PAR(cal, shk)
    return cal, solution_interface(model, PAR, SS_numeric(PAR, tg))
end

irf(model, sol, state, sd) = DataFrame(simulate_model((; model..., ne = 1), sol.sol_mat, T_IR,
                                                      shock_col(state, sd), sol.SS, true, true), varnames)

function cor_uv(model, sol, eta_m)
    sim = DataFrame(simulate_model((; model..., ne = size(eta_m, 2)), sol.sol_mat, T_SM, eta_m,
                                   sol.SS, false, true), varnames)
    q = monthly_to_quarterly(sim[!, [:u, :v]])
    hu, hv = hp_filter(Float64.(q.u), LAM), hp_filter(Float64.(q.v), LAM)
    return cor(hu, hv), 100 * std(hu), 100 * std(hv)
end

# ── Part A ───────────────────────────────────────────────────────────────────
log!("Part A — baseline at the code defaults over Λ̄_r. Responses in % (100 × log dev.), month 1 = impact.")
log!("Simulated moments: $T_SM months, quarterly averages, HP λ = $(Int(LAM)) on log deviations; data cor(u,v) = −0.804.")
model = build_model(TARGETS)
S = SHOCKS_BASE
eta_all = hcat(shock_col(:z, S.σz), shock_col(:δ, S.σδ), shock_col(:s, S.σs))
for Λr in GRID
    tg = merge(TARGETS, (Λr_target = Λr,))
    cal, sol = solve(model, tg)
    tag = isnothing(Λr) ? "costless" : @sprintf("Λ̄_r = %.2f", Λr)
    log!(@sprintf("\n[%s]  α = %.4g", tag, cal.α))
    for (st, sd) in ((:z, S.σz), (:δ, S.σδ), (:s, S.σs))
        r = 100 .* irf(model, sol, st, sd)
        log!(@sprintf("  %s shock: u[1,12] = %+7.3f %+7.3f   v[1,12] = %+7.3f %+7.3f", st,
                      r.u[2], r.u[13], r.v[2], r.v[13]))
    end
    c, su, sv = cor_uv(model, sol, eta_all)
    log!(@sprintf("  all shocks: cor(u,v) = %+.3f   σ(u) = %.2f%%   σ(v) = %.2f%%", c, su, sv))
    for (j, st) in enumerate((:z, :δ, :s))
        cj, _, _ = cor_uv(model, sol, eta_all[:, j:j])
        log!(@sprintf("  %s only:    cor(u,v) = %+.3f", st, cj))
    end
end

write("R8_step7_baseline.txt", join(lines, "\n") * "\n")
