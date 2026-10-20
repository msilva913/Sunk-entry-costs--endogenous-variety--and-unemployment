# Q2_p0_floor.jl
# =============================================================================
# Exploratory, NOT draft-bound. How far must p_0 rise for the reposting margin to deliver
# the data's cor(u, v) = −0.804? Λ_r = F(Q/α) ≥ 1 − p_0, so p_0 sets the floor of Λ̄_r.
#
# Grid: p_0 ∈ {0.5, …, 0.95} × Λ̄_r ∈ {0.55, 0.45, …, 0.10}, kept only when
# Λ̄_r > 1 − p_0 + 0.03. Everything else at TARGETS (PATH B: ψ re-pinned by the exit
# elasticity at each p_0). Shocks: part6c (a), the trivariate VAR, orthogonalized.
# Moments as R8_step7_baseline.jl. Data: cor(u,v) −0.804, σ(u) 14.96 %, σ(v) 14.19 %.
#
# Output: Q2_p0_floor.txt
# =============================================================================

include("run_solution_core.jl")
using Printf, Statistics, CSV

const T_SM  = 60_000
const LAM   = 1_600.0
const P0    = (0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
const LRGRD = (0.55, 0.45, 0.35, 0.25, 0.15, 0.10)

lines = String[]
log!(str) = (push!(lines, str); println(str))

const CAND = CSV.read(joinpath(@__DIR__, "..", "Bartek analysis", "data", "results",
                               "shock_var_trivariate.csv"), DataFrame)
function shocks_from(spec)
    get(sh) = (r = CAND[(CAND.spec .== spec) .& (CAND.shock .== sh), :]; (r.rho_m[1], r.sigma_m[1]))
    (ρz, σz), (ρδ, σδ), (ρs, σs) = get("z"), get("delta"), get("s")
    return (; ρz, σz, ρδ, σδ, ρs, σs)
end
const SH = shocks_from("a_orthogonalized")

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

function moments_uv(sol, eta_m)
    sim = DataFrame(simulate_model((; model..., ne = size(eta_m, 2)), sol.sol_mat, T_SM, eta_m,
                                   sol.SS, false, true), varnames)
    q = monthly_to_quarterly(sim[!, [:u, :v]])
    hu, hv = hp_filter(Float64.(q.u), LAM), hp_filter(Float64.(q.v), LAM)
    return cor(hu, hv), 100 * std(hu), 100 * std(hv)
end

eta_all = hcat(shock_col(:z, SH.σz), shock_col(:δ, SH.σδ), shock_col(:s, SH.σs))
log!("Q2 p_0 floor — cor(u,v) over (p_0, Λ̄_r); shocks part6c (a); $T_SM months, quarterly, HP λ = $(Int(LAM)).")
log!("Data: cor(u,v) −0.804, σ(u) 14.96 %, σ(v) 14.19 %.")
for p0 in P0
    log!(@sprintf("\n== p_0 = %.2f  (floor 1 − p_0 = %.2f)", p0, 1 - p0))
    for Λr in LRGRD
        Λr > 1 - p0 + 0.03 || continue
        tg = merge(TARGETS, (p_0 = p0, Λr_target = Λr))
        res = try
            cal = calibrate_shares(tg)
            PAR = build_PAR(cal, SH)
            sol = solution_interface(model, PAR, SS_numeric(PAR, tg))
            (cal, sol)
        catch err
            log!(@sprintf("  Λ̄_r = %.2f  failed: %s", Λr, first(sprint(showerror, err), 110)))
            nothing
        end
        isnothing(res) && continue
        cal, sol = res
        c, su, sv = moments_uv(sol, eta_all)
        cs = moments_uv(sol, eta_all[:, 3:3])[1]
        log!(@sprintf("  Λ̄_r = %.2f   cor(u,v) %+.3f   s only %+.3f   σ(u) %5.2f%%   σ(v) %5.2f%%   ψ %.4f   log10 α %.0f",
                      Λr, c, cs, su, sv, cal.ψ, log10(cal.α)))
    end
end

write(joinpath(@__DIR__, "Q2_p0_floor.txt"), join(lines, "\n") * "\n")
