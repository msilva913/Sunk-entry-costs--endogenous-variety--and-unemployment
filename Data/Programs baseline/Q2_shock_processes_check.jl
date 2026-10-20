# Q2_shock_processes_check.jl
# =============================================================================
# Q2 model check. Exploratory, NOT draft-bound (§5.3 is regenerated in Q5).
# Simulated cor(u, v) and σ(u) over Λ̄_r at two sets of shock processes:
#   part6b       SHOCKS_BASE: bivariate (z, δ) VAR + univariate AR(1) for s
#   part6c (a)   trivariate (z, δ, s) VAR, Cholesky z → δ → s
#                (Bartek analysis/data/results/shock_var_trivariate.csv)
# Moments as R8_step7_baseline.jl: monthly → quarterly averages → HP λ = 1,600.
# Data: cor(u, v) = −0.804. The grid stops above the floor Λ̄_r > 1 − p_0 = 0.5.
#
# Output: Q2_shock_processes_check.txt
# =============================================================================

include("run_solution_core.jl")
using Printf, Statistics, CSV

const T_SM = 60_000
const LAM  = 1_600.0
const GRID = (nothing, 0.85, 0.70, 0.60, 0.55)

lines = String[]
log!(str) = (push!(lines, str); println(str))

# Candidate processes written by part6c (monthly ρ, σ).
const CAND = CSV.read(joinpath(@__DIR__, "..", "Bartek analysis", "data", "results",
                               "shock_var_trivariate.csv"), DataFrame)
function shocks_from(spec)
    get(sh) = (r = CAND[(CAND.spec .== spec) .& (CAND.shock .== sh), :]; (r.rho_m[1], r.sigma_m[1]))
    (ρz, σz), (ρδ, σδ), (ρs, σs) = get("z"), get("delta"), get("s")
    return (; ρz, σz, ρδ, σδ, ρs, σs)
end
const SETS = (("part6b", SHOCKS_BASE), ("part6c (a)", shocks_from("a_orthogonalized")))

model = let
    m = (parameters = parameters, estimate = estimate, estimation = position,
         npar = length(parameters), ns = length(estimate), priors = priors,
         x = x, y = y, xp = xp, yp = yp, variables = variables, varnames = varnames,
         nx = nx, ny = ny, nvar = nvar, e = ex, eta = eta, ne = ne,
         f = gen_model_equations(), nf = nvar, SS = SS_symbolics(parameters, TARGETS),
         PAR_SS = parameters[:], flag_order = flag_order, flag_deviation = flag_deviation,
         flag_SSsolver = flag_SSsolver)
    process_model(m)
    m
end

function cor_uv(sol, eta_m)
    sim = DataFrame(simulate_model((; model..., ne = size(eta_m, 2)), sol.sol_mat, T_SM, eta_m,
                                   sol.SS, false, true), varnames)
    q = monthly_to_quarterly(sim[!, [:u, :v]])
    hu, hv = hp_filter(Float64.(q.u), LAM), hp_filter(Float64.(q.v), LAM)
    return cor(hu, hv), 100 * std(hu)
end

log!("Q2 model check — simulated cor(u,v) over Λ̄_r; $T_SM months, quarterly averages, HP λ = $(Int(LAM)).")
log!("Data cor(u,v) = −0.804. Columns: all shocks, then each shock alone.")
for (lab, S) in SETS
    log!(@sprintf("\n== %s   ρz %.4f σz %.4f | ρδ %.4f σδ %.4f | ρs %.4f σs %.4f",
                  lab, S.ρz, S.σz, S.ρδ, S.σδ, S.ρs, S.σs))
    eta_all = hcat(shock_col(:z, S.σz), shock_col(:δ, S.σδ), shock_col(:s, S.σs))
    for Λr in GRID
        tg  = merge(TARGETS, (Λr_target = Λr,))
        PAR = build_PAR(calibrate_shares(tg), S)
        sol = solution_interface(model, PAR, SS_numeric(PAR, tg))
        c, su = cor_uv(sol, eta_all)
        cz, cδ, cs = (cor_uv(sol, eta_all[:, j:j])[1] for j in 1:3)
        log!(@sprintf("  Λ̄_r = %-8s cor(u,v) %+.3f   σ(u) %5.2f%%   z only %+.3f   δ only %+.3f   s only %+.3f",
                      isnothing(Λr) ? "costless" : string(Λr), c, su, cz, cδ, cs))
    end
end

write(joinpath(@__DIR__, "Q2_shock_processes_check.txt"), join(lines, "\n") * "\n")
