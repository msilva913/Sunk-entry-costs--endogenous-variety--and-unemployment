# sigma_determinacy.jl
# =============================================================================
# Determinacy of the full model as risk aversion σ falls (Oct 15, 2026; Q1 step 7 follow-up).
#
# Finding that motivated it: at σ ≤ 1e-3 the model has no stable solution, before and after
# R8 (R8_test_free_entry.jl). Question: is the failure the variety externality ρ(N) (strategic
# complementarity in entry), or the business-formation block under near-linear utility?
#
# Configurations (all at the canonical TARGETS: D1, PATH B, costless reposting):
#   A  baseline (p_0 = 0.5, variety on)
#   B  p_0 = 0 (no endogenous exit), variety on
#   C  p_0 = 0, ρ ≡ 1 (f[13] overridden: no variety effect, ζ = 0)
#   D  as C, with the business-formation block frozen: N ≡ 1 (f[22]) and entry labor at its
#      steady-state value (f[4] replaces the firm Euler). This approximates AGS, where
#      clone replacement removes entry as a decision.
# For each σ: the number of stable generalized eigenvalues nk (|T_ii| < |S_ii|, as in
# first_order) against nx = 8, the BK verdict, and the eigenvalue moduli nearest 1.
#
# Output: sigma_determinacy.txt
# =============================================================================

include("run_solution_core.jl")
using Printf

const SHK = (ρz = 0.902, σz = 0.0092, ρδ = 0.592, σδ = 0.0669, ρs = 0.874, σs = 0.0854)
const SIGMAS = (1.0, 0.5, 0.2, 0.1, 0.05, 0.01, 1e-3, 0.0)

tA = TARGETS
tB = merge(TARGETS, (p_0 = 0.0, dest_end_frac = 0.0, Xc_Y = 0.0))

function build_model(tg; novar = false, freeze = false)
    f = gen_model_equations()
    novar  && (f[13] = ρ - 1)
    if freeze
        st = steady_state(calibrate_shares(tg))
        f[22] = N - 1                 # N ≡ 1: no business-formation dynamics
        f[4]  = L_e - st.L_e          # entry labor at its steady state (replaces the firm Euler)
    end
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

# Generalized eigenvalues of the first-order pencil, ordered as in first_order
function bk_report(model, tg)
    cal = calibrate_shares(tg)
    PAR = build_PAR(cal, SHK)
    SS  = SS_numeric(PAR, tg)
    PAR_SS = Base.invokelatest(getfield(Main, :eval_PAR_SS), PAR)
    d   = Base.invokelatest(getfield(Main, :eval_deriv), PAR_SS, SS)
    fd  = d.d1
    fx, fy   = fd[:, 1:nx], fd[:, nx+1:nvar]
    fxp, fyp = fd[:, nvar+1:nvar+nx], fd[:, nvar+nx+1:2nvar]
    F   = schur([-fxp -fyp], [fx fy])
    t, s = abs.(diag(F.T)), abs.(diag(F.S))
    nk  = count(t .< s)
    mods = [s[i] == 0 ? Inf : t[i] / s[i] for i in eachindex(t)]   # |eigenvalue|; < 1 stable
    near = sort(filter(m -> 0.5 < m < 2.0, mods))
    return (; nk, near, ninf = count(isinf, mods))
end

lines = String[]
log!(s) = (push!(lines, s); println(s))
log!("Determinacy as σ falls — stable roots nk (need nx = $nx); |eigenvalues| in (0.5, 2)")

configs = (("A baseline (p_0 = 0.5, variety)",       tA, false, false),
           ("B p_0 = 0, variety",                    tB, false, false),
           ("C p_0 = 0, ρ ≡ 1",                      tB, true,  false),
           ("D p_0 = 0, ρ ≡ 1, business block frozen", tB, true,  true))
for (lab, tg, nv, fr) in configs
    model = build_model(tg; novar = nv, freeze = fr)
    log!("\n== $lab")
    for σv in SIGMAS
        r = try
            bk_report(model, merge(tg, (σ = σv,)))
        catch e
            (; nk = -1, near = Float64[], ninf = -1, err = first(sprint(showerror, e), 90))
        end
        verdict = r.nk == nx ? "determinate" : r.nk > nx ? "INDETERMINATE (too many stable)" :
                  r.nk < 0 ? "ERROR" : "NO STABLE SOLUTION (too few stable)"
        log!(@sprintf("  σ = %-6g nk = %2d  %-34s near-unit |λ|: %s", σv, r.nk, verdict,
                      join([@sprintf("%.4f", m) for m in r.near], " ")))
        hasproperty(r, :err) && log!("      " * r.err)
    end
end
write("sigma_determinacy.txt", join(lines, "\n") * "\n")

# ── Critical σ* by bisection: the smallest σ with a unique bounded solution ──────────────
# Assumes determinacy at σ = 0.5 and failure at σ = 0 (true in A–C above); 20 halvings.
log!("\nCritical σ* (bisection; unique bounded solution for σ > σ*):")
for (lab, tg, nv, fr) in configs[1:3]
    model = build_model(tg; novar = nv, freeze = fr)
    ok(σv) = bk_report(model, merge(tg, (σ = σv,))).nk == nx
    lo, hi = 0.0, 0.5
    (ok(hi) && !ok(lo)) || (log!("  $lab: bracket fails"); continue)
    for _ in 1:20
        mid = (lo + hi) / 2
        ok(mid) ? (hi = mid) : (lo = mid)
    end
    log!(@sprintf("  %-40s σ* ∈ (%.5f, %.5f)", lab, lo, hi))
end
write("sigma_determinacy.txt", join(lines, "\n") * "\n")
