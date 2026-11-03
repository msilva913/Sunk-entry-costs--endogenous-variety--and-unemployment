# Q2_amplification_sweep.jl
# =============================================================================
# Exploratory, NOT draft-bound, NOT a calibration. Which existing margins amplify the z and
# δ shocks? One parameter at a time around a working point; nothing is chosen here (values
# come from Q3/Q6). η_L and ε are data-pinned and are not swept (MS, Oct 20).
#
# Levers: b/w; x_v (sunk vs per-hire cost κ, Pissarides 2009); X_Y (recruiting cost share,
# which together with b/w and x_v pins the bargaining power ϕ in Stage 4); ξ_inv (vacancy
# creation elasticity); δ-specific: dest_elast_target (exit elasticity), dest_end_frac (ω_δ).
#
# Working point: TARGETS + p_0 = 0.8, Λ̄_r = 0.25 (low reposting, as the data SVAR's s → v
# response suggests; Q2_p0_floor.txt). Shocks: SHOCKS_BASE (S14 values since Oct 29; part6c (a)
# before). A costless-reposting row is shown for reference.
#
# Fundamental surplus (after Ljungqvist and Sargent 2017), steady state, per period, as a
# share of the recruiter's revenue per worker w_int:
#   FS0 = (w_int − b)/w_int
#   FS  = (w_int − b − (r+τ)κ − s𝓡)/w_int,   τ = δ_e + (1−δ_e)s,  𝓡 = Q − K − H
# (r+τ)κ annuitizes the per-hire cost over the expected match duration; s𝓡 is the reposting
# shortfall that enters the Nash wage. Amplification of z scales roughly with 1/FS.
#
# Reported: ϕ, FS0, FS; cor(u,v), σ(u), σ(v), σ(θ)/σ(LP) with all shocks; σ(u) and the share
# of var(u) for each shock alone (independent shocks: shares add up to about 1).
# Data, 1951–2019 (S13): cor −0.925, σ(u) 12.48 %, σ(v) 13.82 % (findings.md "Sample windows").
# σ(θ)/σ(LP) 11.70 is the old 1951–2024 value. Data SVAR shares of var(u) (scratch, Oct 20,
# not yet a program): s ≤ ~25 %, z + δ ≈ 55–60 % with shocks ordered first.
#
# δ-process block (Nov 2, 2026): an upper-bound test at the working point. How much can a more
# persistent δ move u? Cases: the HP-matched δ (SHOCKS_BASE); the Hamilton-matched δ from part6b
# (shock_calibration.csv, hamilton_matched_*); and ρ_δ raised with σ_δ rescaled so δ's unconditional
# variance stays at the baseline value (persistence alone). Exploratory, not a calibration.
#
# Clone-replacement block (Nov 3, 2026): the working point with only the product-line entry Euler
# f[4] replaced by Ñ_{t+1} = 1 (entry replaces exactly what δ destroys, unpriced), as in R8 Test 6.
# Everything else, including the costly reposting margin (Λ̄_r = 0.25), endogenous exit (p_0 = 0.8)
# and variety, is unchanged. It isolates the replacement-pricing channel (R19).
#
# Output: Q2_amplification_sweep.txt
# =============================================================================

include("run_solution_core.jl")
using Printf, Statistics, CSV

const T_SM  = 60_000
const LAM   = 1_600.0
const BASE  = merge(TARGETS, (p_0 = 0.8, Λr_target = 0.25))
const SWEEP = ((:b_ratio, (0.8, 0.95)), (:x_v, (0.2, 0.9)), (:X_Y, (0.0075, 0.03)),
               (:ξ_inv, (0.25, 4.0)), (:dest_elast_target, (2.5, 10.0)),
               (:dest_end_frac, (0.25, 0.75)))

lines = String[]
log!(str) = (push!(lines, str); println(str))

const SH  = SHOCKS_BASE
const ETA = hcat(shock_col(:z, SH.σz), shock_col(:δ, SH.σδ), shock_col(:s, SH.σs))

function build_model(; clone_N = nothing)
    f = gen_model_equations()
    # Clone replacement: the product-line entry Euler gives way to Ñ_{t+1} = Ñ̄, so entry replaces
    # exactly what δ destroys. Ñ̄ is the working point's steady-state pre-exit stock (AGS's
    # normalization Ñ = 1 does not hold here).
    clone_N === nothing || (f[4] = Ntp - clone_N)
    m = (parameters = parameters, estimate = estimate, estimation = position,
         npar = length(parameters), ns = length(estimate), priors = priors,
         x = x, y = y, xp = xp, yp = yp, variables = variables, varnames = varnames,
         nx = nx, ny = ny, nvar = nvar, e = ex, eta = eta, ne = ne,
         f = f, nf = nvar, SS = SS_symbolics(parameters, TARGETS),
         PAR_SS = parameters[:], flag_order = flag_order, flag_deviation = flag_deviation,
         flag_SSsolver = flag_SSsolver)
    process_model(m)
    m
end
model = build_model()

hpq(sim, v) = hp_filter(Float64.(monthly_to_quarterly(sim[!, [v]])[!, v]), LAM)

function stats(sol, eta_m; mdl = model)
    sim = DataFrame(simulate_model((; mdl..., ne = size(eta_m, 2)), sol.sol_mat, T_SM, eta_m,
                                   sol.SS, false, true), varnames)
    hu, hv, hθ, hlp = hpq(sim, :u), hpq(sim, :v), hpq(sim, :θ), hpq(sim, :labor_prod)
    return (c = cor(hu, hv), su = 100std(hu), sv = 100std(hv), amp = std(hθ) / std(hlp))
end

function fundamental_surplus(cal, ss)
    τ  = ss.δ_e + (1 - ss.δ_e) * cal.s
    sh = ss.Q - ss.K - ss.H                                   # reposting shortfall 𝓡
    return (ss.w_int - cal.b) / ss.w_int,
           (ss.w_int - cal.b - (cal.r + τ) * cal.κ - cal.s * sh) / ss.w_int
end

function run_case(label, tg; shk = SH, mdl = model)
    eta = hcat(shock_col(:z, shk.σz), shock_col(:δ, shk.σδ), shock_col(:s, shk.σs))
    res = try
        cal = calibrate_shares(tg)
        PAR = build_PAR(cal, shk)
        (cal, solution_interface(mdl, PAR, SS_numeric(PAR, tg)))
    catch err
        log!(@sprintf("  %-26s failed: %s", label, first(sprint(showerror, err), 100)))
        return
    end
    cal, sol = res
    fs0, fs = fundamental_surplus(cal, sol.ss)
    a   = stats(sol, eta; mdl)
    one = [stats(sol, eta[:, j:j]; mdl) for j in 1:3]
    v2  = [o.su^2 for o in one]
    sh  = 100 .* v2 ./ sum(v2)
    log!(@sprintf("  %-26s ϕ %.3f  FS0 %.3f  FS %.3f | cor %+.3f  σu %5.2f  σv %5.2f  σθ/σLP %5.2f | σu by z %5.2f δ %5.2f s %5.2f | share z %4.1f δ %4.1f s %4.1f",
                  label, cal.ϕ, fs0, fs, a.c, a.su, a.sv, a.amp, one[1].su, one[2].su, one[3].su, sh...))
end

log!(@sprintf("Q2 amplification sweep — working point p_0 = 0.8, Λ̄_r = 0.25; shocks SHOCKS_BASE (S14): ρz %.4f σz %.4f ρδ %.4f σδ %.4f ρs %.4f σs %.4f; σ in %%, HP-1600 quarterly.", SH...))
log!("Data 1951–2019: cor −0.925, σu 12.48, σv 13.82 (σθ/σLP 11.70 is 1951–2024); SVAR shares of var(u) (scratch): s ≤ ~25 %, z+δ ≈ 55–60 %.")
log!("\n== reference")
run_case("costless reposting", merge(TARGETS, (p_0 = 0.8,)))
run_case("working point", BASE)
for (par, vals) in SWEEP
    log!(@sprintf("\n== %s  (working point %s)", par, string(getfield(BASE, par))))
    for val in vals
        run_case(@sprintf("%s = %s", par, string(val)), merge(BASE, NamedTuple{(par,)}((val,))))
    end
end

# ── δ process: upper-bound test at the working point ───────────────────────
const CAL = CSV.read(joinpath(@__DIR__, "..", "Bartek analysis", "data", "results",
                              "shock_calibration.csv"), DataFrame)
const DROW = CAL[CAL.series .== "delta", :]
const ρδ_ham, σδ_ham = DROW.hamilton_matched_rho_m[1], DROW.hamilton_matched_sigma_m[1]
same_var(ρ) = SH.σδ * sqrt((1 - ρ^2) / (1 - SH.ρδ^2))   # keeps σδ²/(1-ρδ²) at the baseline
log!(@sprintf("\n== δ process (upper-bound test; working point; baseline ρδ %.4f σδ %.4f)", SH.ρδ, SH.σδ))
for (label, ρ, σ) in (("δ HP-matched (baseline)", SH.ρδ, SH.σδ),
                      (@sprintf("δ Hamilton-matched %.3f", ρδ_ham), ρδ_ham, σδ_ham),
                      (@sprintf("δ ρ %.3f, base variance", ρδ_ham), ρδ_ham, same_var(ρδ_ham)),
                      ("δ ρ 0.970, base variance", 0.97, same_var(0.97)))
    run_case(@sprintf("%s (σδ %.4f)", label, σ), BASE; shk = merge(SH, (ρδ = ρ, σδ = σ)))
end

# ── clone replacement at the working point (R19 isolated) ──────────────────
# process_model compiles the equations into global functions, so building a second model
# replaces the first: this block must come last, and its reference row is "working point" above.
function pre_exit_N(tg)
    cal = calibrate_shares(tg); st = steady_state(cal)
    return st.N / ((1 - st.δ_e) / (1 - cal.δ))
end
const NT_SS = pre_exit_N(BASE)
const MODEL_CLONE = build_model(clone_N = NT_SS)
log!(@sprintf("\n== clone replacement (f[4]: Ñ_{t+1} = Ñ̄ = %.6f; working point, matched shocks; compare \"working point\" above)", NT_SS))
run_case("clone replacement", BASE; mdl = MODEL_CLONE)
run_case("clone repl., costless repost", merge(TARGETS, (p_0 = 0.8,)); mdl = MODEL_CLONE)

write(joinpath(@__DIR__, "Q2_amplification_sweep.txt"), join(lines, "\n") * "\n")
