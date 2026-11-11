# Q2_persistence_comparison.jl
# =============================================================================
# Exploratory, not draft-bound. TRUE-SHOCK statistics (decisions.md S16): models against each
# other, not against the data (that needs the full model through the SVAR, gate 2).
#
# Question: does the full model's business-formation block add persistence or Beveridge-curve
# comovement relative to AGS (Programs AGS/ags_persistence.jl)? The SVAR comparison found that in
# the data every shock moves u persistently (ρ₁ ≈ 0.95) and against v (cor ≈ −0.97), while AGS
# does so only for z (q3_model_svar.txt §5).
#
# Full-model cases, all with AGS's reactivation (Λ̄_r 0.25, p_r 0.8, ψ_r 1; decoupled from exit)
# and the shared targets (x_v 0.5, ξ_inv 1, b/w 0.9, X/Y 1.5 %), with Block B′'s shock processes:
#   no exit      p_0 = 0 (δ = measured exit, as in AGS): variety and priced product-line entry
#                are the only additions to AGS (and σ = 1 against AGS's 0: determinacy needs
#                σ > 0.011 under priced entry, sigma_determinacy.txt);
#   working pt   p_0 = 0.8, half of measured exit endogenous (TARGETS): exogenous δ̄ is half of
#                AGS's, so the same proportional δ shock is half as large in absolute terms;
# each with priced entry (the entry Euler) and clone replacement (f[4]: Ñ_{t+1} = Ñ̄; isolates
# the pricing of replacement, R19).
# σ check (Nov 6): without exit and with clone replacement, the full model's δ effect on u is twice
# AGS's although both hold product lines fixed. Remaining differences: σ = 1 (δ moves consumption
# and the discount factor) and labor used to replace destroyed lines. The no-exit cases are rerun
# at σ = 0.1 and 0.03 (clone) and 0.1 (priced); what is left at small σ is the entry-labor channel.
# For each shock alone: σ(u), σ(v), ρ₁(u), ρ₁(v), cor(u, v), cor(u, labor productivity) on
# quarterly averages, HP-1600, and the quarter in which u's impulse response peaks. Statistics
# from time_series_fun.jl (shock_moments, peak_quarter), as for AGS.
#
# Output: Q2_persistence_comparison.txt. Run Programs AGS/ags_persistence.jl first.
# =============================================================================

include("run_solution_core.jl")
using Printf, Statistics, CSV

const SHOCK_CSV = joinpath(@__DIR__, "..", "Bartek analysis", "data", "results", "shock_calibration.csv")
const AGS_CSV   = joinpath(@__DIR__, "..", "Programs AGS", "ags_persistence.csv")
const SHK_B     = shock_processes(SHOCK_CSV; prefix = "bprime")
const T_SM, LAM, T_IRF_M = 60_000, 1_600.0, 60
const REPOST    = (Λr_target = 0.25, p_r = 0.8, ψ_r = 1.0)
const NO_EXIT = merge(TARGETS, (p_0 = 0.0, dest_end_frac = 0.0, Xc_Y = 0.0), REPOST)
const CASES = (
    ("no exit",    NO_EXIT),
    ("working pt", merge(TARGETS, (p_0 = 0.8,), REPOST)),
)
const SIGMA_PRICED = (("no exit, σ 0.1", merge(NO_EXIT, (σ = 0.1,))),)
const SIGMA_CLONE  = (("no exit, σ 0.1", merge(NO_EXIT, (σ = 0.1,))),
                      ("no exit, σ 0.03", merge(NO_EXIT, (σ = 0.03,))))

lines = String[]
log!(str) = (push!(lines, str); println(str))

function build_model(; clone_N = nothing)
    f = gen_model_equations()
    if !isnothing(clone_N)
        f[4] = Ntp - clone_N            # clone replacement: entry replaces exactly what δ destroys
    end
    m = (parameters = parameters, estimate = estimate, estimation = position,
         npar = length(parameters), ns = length(estimate), priors = priors,
         x = x, y = y, xp = xp, yp = yp, variables = variables, varnames = varnames,
         nx = nx, ny = ny, nvar = nvar, e = ex, eta = eta, ne = ne,
         f = f, nf = nvar, SS = SS_symbolics(parameters, TARGETS),
         PAR_SS = parameters[:], flag_order = flag_order, flag_deviation = flag_deviation,
         flag_SSsolver = flag_SSsolver)
    process_model(m)
    return m
end

"""Pre-exit steady-state stock of product lines Ñ̄ (clone replacement holds it fixed)."""
function pre_exit_N(tg)
    cal = calibrate_shares(tg)
    st = steady_state(cal)
    return st.N / ((1 - st.δ_e) / (1 - cal.δ))
end

"""Quarterly HP cycles of u, v and measured labor productivity from a monthly simulation."""
function cycles(sim)
    q = monthly_to_quarterly(DataFrame(u = sim.u, v = sim.v, lp = sim.labor_prod))
    return DataFrame([c => hp_filter(Float64.(q[!, c]), LAM) for c in names(q)])
end

function run_case(label, tg, mdl)
    cal = calibrate_shares(tg)
    PAR = build_PAR(cal, SHK_B)
    sol = solution_interface(mdl, PAR, SS_numeric(PAR, tg))
    out = DataFrame()
    for (state, name, sd) in ((:z, "z", SHK_B.σz), (:δ, "δ", SHK_B.σδ), (:s, "s", SHK_B.σs))
        col = shock_col(state, sd)
        sim = DataFrame(simulate_model((; mdl..., ne = 1), sol.sol_mat, T_SM, col, sol.SS,
                                       false, true), varnames)
        mom = shock_moments(cycles(sim); prod = "lp")
        ir = DataFrame(simulate_model((; mdl..., ne = 1), sol.sol_mat, T_IRF_M + 1, col, sol.SS,
                                      true, true), varnames)
        kq, val = peak_quarter(100 .* ir.u[2:end])
        push!(out, (; case = label, shock = name, mom..., u_peak_q = kq, u_peak = val))
    end
    return out
end

# Priced entry first, then clone replacement: process_model redefines the global evaluation
# functions, so each model is built just before its cases run.
res = DataFrame()
model = build_model()
for (lab, tg) in (CASES..., SIGMA_PRICED...)
    append!(res, run_case("full, $lab, priced", tg, model))
end
for (lab, tg) in (CASES..., SIGMA_CLONE...)
    m_clone = build_model(clone_N = pre_exit_N(tg))
    append!(res, run_case("full, $lab, clone", tg, m_clone))
end
append!(res, CSV.read(AGS_CSV, DataFrame); cols = :union)

log!("Persistence and comovement by shock: AGS against the full model (TRUE shocks, not comparable with")
log!("the data's SVAR numbers; decisions.md S16). Block B′ shock processes: " *
     @sprintf("ρz %.4f σz %.4f ρδ %.4f σδ %.4f ρs %.4f σs %.4f (monthly).", SHK_B...))
log!("Reactivation as in AGS (Λ̄_r 0.25, p_r 0.8, ψ_r 1). Quarterly averages, HP-1600; u peak: quarter of the")
log!("largest |response| of u to a one-SD shock (0 = impact), and its value in %.")
log!("Data, through SVAR A (q3_model_svar.txt §5): ρ₁(u) 0.96 / 0.96 / 0.94 and cor(u,v) −0.98 / −0.98 / −0.96")
log!("for z / δ / s; u peaks at h 4 / 4 / 2.")
order = ["AGS", "full, no exit, priced", "full, no exit, σ 0.1, priced", "full, no exit, clone",
         "full, no exit, σ 0.1, clone", "full, no exit, σ 0.03, clone", "full, working pt, priced",
         "full, working pt, clone"]
for name in ("z", "δ", "s")
    log!("\n$name shock alone")
    log!(@sprintf("  %-30s %6s %6s %6s %6s %7s %8s  %s", "model", "σ(u)", "σ(v)", "ρ₁(u)", "ρ₁(v)",
                  "cor(u,v)", "cor(u,LP)", "u peak (quarter)"))
    for c in order
        r = res[(res.case .== c) .& (res.shock .== name), :]
        if nrow(r) == 0
            continue
        end
        r = r[1, :]
        log!(@sprintf("  %-30s %6.2f %6.2f %6.3f %6.3f %+7.3f %+8.3f  %+.3f (%d)", c, r.σu, r.σv, r.ρu,
                      r.ρv, r.cuv, r.cup, r.u_peak, r.u_peak_q))
    end
end
log!("\ncor(u,LP): AGS's labor productivity is z; the full model's is measured Y/(ρL). NaN: LP does not move.")

open("Q2_persistence_comparison.txt", "w") do io
    foreach(l -> println(io, l), lines)
end
