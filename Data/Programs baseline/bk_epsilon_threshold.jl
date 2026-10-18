# bk_epsilon_threshold.jl
# =============================================================================
# Blanchard–Kahn threshold in ε, and which feedback loop sets it.
#
# Hypothesis (Oct 18, 2026). Under DS-CES, goods profit per line ∝ N^{(2−ε)/(ε−1)}:
# dilution (1/N) against the variety effect (ρ(N) = N^{1/(ε−1)}). Entry costs are paid in
# effective labor, so ν_f = ρ f_e/μ scales with ρ and the variety effect cancels on the
# entry margin. Costs paid in goods (continuation draws χ, vacancy costs G(e) and κ) do
# not scale with ρ, which opens two loops of gain ∝ 1/(ε−1):
#   exit loop:          N ↓ → profit per line ↓ (ε < 2) → χ^c ↓ → more exit → N ↓
#   labor-market loop:  N ↓ → w_int = ρ(N) z/μ ↓ → vacancies ↓ → L^c ↓ → profit ↓ → entry ↓
# A predetermined root then crosses the unit circle and no bounded equilibrium exists
# (qzflag = 0: more explosive roots than jumps). The threshold ε* is found by bisection
# under four specifications:
#   base     code defaults (D1, PATH B, p_0 = 0.5, b/w = 0.9, x_v = 0.5, costless reposting)
#   no_exit  p_0 = 0 (dest_end_frac = 0, Xc_Y = 0): closes the exit loop
#   low_b    b/w = 0.71: larger fundamental surplus, weakens the labor-market loop
#   no_var   ρ ≡ 1 (f[13] replaced, as in Comparison C): removes the variety effect.
#            With N̄ = 1, ρ(N̄) = 1, so the steady state is unchanged.
#
# Output: bk_epsilon_threshold.txt
# =============================================================================

include("run_solution_core.jl")
using Printf, LinearAlgebra

const SPECS = [
    (:base,    NamedTuple(),                                     false),
    (:no_exit, (p_0 = 0.0, dest_end_frac = 0.0, Xc_Y = 0.0),     false),
    (:low_b,   (b_ratio = 0.71,),                                false),
    (:no_var,  NamedTuple(),                                     true),
]

lines = String[]
log!(s) = (push!(lines, s); println(s))

function build_model(tg, no_var)
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

# (determinate?, dominant |eigenvalue| of hx or NaN)
function check(model, tg)
    cal = calibrate_shares(tg)
    PAR = build_PAR(cal, SHOCKS_BASE)
    try
        sol = solution_interface(model, PAR, SS_numeric(PAR, tg))
        return true, maximum(abs.(eigvals(sol.sol_mat.hx)))
    catch err
        occursin("Blanchard", sprint(showerror, err)) || rethrow()
        return false, NaN
    end
end

log!("BK threshold in ε by specification (bisection to 0.005)")
results = Dict{Symbol,Any}()
for (name, over, no_var) in SPECS
    tg0 = merge(TARGETS, over)
    model = build_model(tg0, no_var)
    grid = Dict{Float64,Tuple{Bool,Float64}}()
    chk(ε) = get!(grid, ε) do; check(model, merge(tg0, (ε = ε,))); end
    log!("\n[$name]")
    for ε in (4.3, 3.0, 2.0, 1.5, 1.2)
        d, ev = chk(ε)
        log!(@sprintf("  ε = %.3f: %s   max|eig(hx)| = %s", ε, d ? "determinate" : "NO bounded equilibrium",
                      isnan(ev) ? "—" : @sprintf("%.6f", ev)))
    end
    pts = sort(collect(keys(grid)))
    fails = [ε for ε in pts if !grid[ε][1]]
    oks   = [ε for ε in pts if grid[ε][1]]
    if isempty(fails)
        results[name] = (εstar = NaN, note = "determinate at every ε tested (down to 1.2)")
    elseif isempty(oks) || minimum(oks) < maximum(fails)
        results[name] = (εstar = NaN, note = "non-monotone determinacy region; inspect grid")
    else
        lo, hi = maximum(fails), minimum(oks)
        while hi - lo > 0.005
            mid = round((lo + hi) / 2; digits = 4)
            chk(mid)[1] ? (hi = mid) : (lo = mid)
        end
        results[name] = (εstar = (lo + hi) / 2, note = @sprintf("fails at %.4f, determinate at %.4f (max|eig| %.6f)",
                                                              lo, hi, grid[hi][2]))
    end
    log!(@sprintf("  ⇒ ε* ≈ %s  (%s)", isnan(results[name].εstar) ? "—" : @sprintf("%.3f", results[name].εstar),
                  results[name].note))
end

log!("\nSummary: ε* by specification")
for (name, _, _) in SPECS
    r = results[name]
    log!(@sprintf("  %-8s ε* = %-7s %s", name, isnan(r.εstar) ? "—" : @sprintf("%.3f", r.εstar), r.note))
end
write("bk_epsilon_threshold.txt", join(lines, "\n") * "\n")
