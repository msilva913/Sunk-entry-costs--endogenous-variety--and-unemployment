# R8_test_reposting.jl
# =============================================================================
# Q1 step 4 tests (Notes/Q1_R8_D1_PATHB_plan.md): reposting dynamics in
# run_solution_core.jl (Λ_r, Q_rep, H; 𝓡 = Q − K − H in the JCC and the wage; X^r).
#
# Test 4 (exact). Costless reposting (no Λr_target ⇒ rp = 0) reproduces the step-2
#   IRFs (R8_test_timing.jls) for both cases.
# Test 5. Partial reposting (Λ̄_r ∈ {0.95, 0.85, 0.70}, base case): the model solves,
#   no unstable roots, and positions are conserved up to exit AND the retired mass:
#   (v_pre + L)_{t+1} = Λ_{t+1}(1−δ_t)[(v_t + L_t) − (1 − Λ_r,t+1) s_t (1 − u_t)],
#   checked at 1e-3 SD, where a correct model leaves only second-order error.
# Report. Separation-shock responses of u and v by Λ̄_r (not draft-bound; Q5).
#
# Output: R8_test_reposting.txt, R8_test_reposting.jls
# =============================================================================

include("run_solution_core.jl")
using Printf

const SHK = (ρz = 0.902, σz = 0.0092, ρδ = 0.592, σδ = 0.0669, ρs = 0.874, σs = 0.0854)
const T_IR = 60

base_targets = merge(TARGETS, (dest_ann = 0.0320, b_ratio = 0.9, x_v = 0.5))
tg_p0  = merge(base_targets, (p_0 = 0.0, dest_end_frac = 0.0, Xc_Y = 0.0))
tg_bas = merge(base_targets, (p_0 = 0.5, dest_end_frac = 0.5, dest_elast_target = 5.0))

function solve_case(tg)
    cal = calibrate_shares(tg)
    PAR = build_PAR(cal, SHK)
    model = (parameters = parameters, estimate = estimate, estimation = position,
             npar = length(parameters), ns = length(estimate), priors = priors,
             x = x, y = y, xp = xp, yp = yp, variables = variables, varnames = varnames,
             nx = nx, ny = ny, nvar = nvar, e = ex, eta = eta, ne = ne,
             f = gen_model_equations(), nf = nvar,
             SS = SS_symbolics(parameters, tg), PAR_SS = parameters[:],
             flag_order = flag_order, flag_deviation = flag_deviation,
             flag_SSsolver = flag_SSsolver)
    process_model(model)
    sol = solution_interface(model, PAR, SS_numeric(PAR, tg))
    m1 = (; model..., ne = 1)
    function irf(state, sd)
        col = zeros(nx, 1); col[findfirst(==(state), Symbol.(x)), 1] = sd
        100 .* DataFrame(simulate_model(m1, sol.sol_mat, T_IR, col, sol.SS, true, true), varnames)
    end
    nunstable = count(abs.(eigvals(sol.sol_mat.hx)) .>= 1 - 1e-10)
    return (; ss = sol.ss, cal, nunstable,
            irf_z = irf(:z, SHK.σz), irf_δ = irf(:δ, SHK.σδ), irf_s = irf(:s, SHK.σs))
end

lines = String[]
log!(s) = (push!(lines, s); println(s))
log!("R8 step 4 tests — reposting dynamics. IRFs in 100 × log deviations, T = $T_IR")

# ── Test 4 ───────────────────────────────────────────────────────────────────
step2 = deserialize("R8_test_timing.jls")
worst4 = 0.0
for (k, tg) in ((:p0zero, tg_p0), (:base, tg_bas))
    r = solve_case(tg)
    for shk in (:irf_z, :irf_δ, :irf_s)
        o, n = getfield(step2[k], shk), getfield(r, shk)
        vars = intersect(Symbol.(names(o)), Symbol.(names(n)))
        d = maximum(maximum(abs.(n[!, v] .- o[!, v])) for v in vars)
        global worst4 = max(worst4, d)
        log!(@sprintf("  Test 4  %-6s %-6s max |costless − step 2| = %.2e", k, shk, d))
    end
end
log!(@sprintf("Test 4 %s (worst %.2e, tolerance 1e-10)", worst4 < 1e-10 ? "PASSED" : "FAILED", worst4))

# ── Test 5 and report ────────────────────────────────────────────────────────
function retire_gap(r, df)
    lv(v) = getfield(r.ss, v) .* exp.(df[!, v] ./ 100)
    v, L, δe, u, Λr = lv(:v), lv(:L), lv(:δ_e), lv(:u), lv(:Λ_r)
    vp   = lv(:v_pret)
    δt   = r.cal.δ .* exp.(df.δ ./ 100)
    δlag = r.cal.δ .* exp.(df.δl ./ 100)
    st   = r.cal.s .* exp.(df.s ./ 100)
    Λ    = (1 .- δe) ./ (1 .- δlag)
    lhs  = vp[2:end] .+ L[2:end]
    rhs  = Λ[2:end] .* (1 .- δt[1:end-1]) .*
           ((v[1:end-1] .+ L[1:end-1]) .- (1 .- Λr[2:end]) .* st[1:end-1] .* (1 .- u[1:end-1]))
    return maximum(abs.(lhs .- rhs))
end

results = Dict{Float64,Any}()
worst5 = 0.0
for Λr in (0.95, 0.85, 0.70)
    r = solve_case(merge(tg_bas, (Λr_target = Λr,)))
    results[Λr] = r
    log!(@sprintf("\nΛ̄_r = %.2f: unstable roots %d (want 0); α = %.4g", Λr, r.nunstable, r.cal.α))
    for shk in (:irf_z, :irf_δ, :irf_s)
        df = getfield(r, shk)
        g1, gs = retire_gap(r, df), retire_gap(r, df ./ 1000)
        global worst5 = max(worst5, gs / g1)
        log!(@sprintf("    position identity with retirement, %-6s gap: 1 SD %.2e, 1e-3 SD %.2e (ratio %.1e)",
                      shk, g1, gs, gs / g1))
    end
    irs = r.irf_s   # not `s`: in the REPL a loop assignment would overwrite the symbolic state s
    log!(@sprintf("    s shock: u[3,13] = %7.4f %7.4f   v[3,13] = %7.4f %7.4f   Λ_r[3] = %8.5f",
                  irs.u[3], irs.u[13], irs.v[3], irs.v[13], irs.Λ_r[3]))
end
log!(@sprintf("\nTest 5 %s (gap ratio at 1e-3 SD: worst %.1e; second order ⇒ ~1e-6)",
              worst5 < 2e-6 ? "PASSED" : "FAILED", worst5))
serialize("R8_test_reposting.jls", results)
write("R8_test_reposting.txt", join(lines, "\n") * "\n")
