# R8_test_timing.jl
# =============================================================================
# Q1 step 2 tests (Notes/Q1_R8_D1_PATHB_plan.md): the timing restructure of
# run_solution_core.jl, with costless reposting.
#
# Test 1 (exact). At p_0 = 0, Λ ≡ 1 and the old and new exit dating coincide, so
#   every variable common to both codes must reproduce R8_snapshot_pre.jls, except
#   δ_e: the old code built δ_e,t from δ_t (the draft's δ_e,t+1), the new code from
#   δ_{t-1}. So new δ_e[t+1] = old δ_e[t].
# Test 2. At p_0 > 0 (case :base): the model solves (SS residual, Blanchard–Kahn),
#   and positions are conserved up to exit in levels:
#   (v_pre + L)_{t+1} = Λ_{t+1}(1-δ_t)(v_t + L_t), Λ_t = (1-δ_e,t)/(1-δ_{t-1}).
#   It also reports how far the base IRFs moved from the pre-R8 code (they should:
#   exit is now dated t+1 in the laws of motion).
#
# Output: R8_test_timing.txt
# =============================================================================

include("run_solution_core.jl")
using Printf

const SHK = (ρz = 0.902, σz = 0.0092, ρδ = 0.592, σδ = 0.0669, ρs = 0.874, σs = 0.0854)
const T_IR = 60

base_targets = merge(TARGETS, (dest_ann = 0.0320, b_ratio = 0.9, x_v = 0.5))
cases = Dict(
    :p0zero => merge(base_targets, (p_0 = 0.0, dest_end_frac = 0.0, Xc_Y = 0.0)),
    :base   => merge(base_targets, (p_0 = 0.5, dest_end_frac = 0.5, dest_elast_target = 5.0)),
)

function solve_case(tg)
    cal = calibrate_shares(tg)
    @unpack f_e, δ, s, z, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, x_m, ψ, f_m, p_0 = cal
    PAR = [f_e; z; δ; s; b; ϕ; r; σ; ε; A; η_L; κ; ξ_inv; x_m; ψ; f_m; p_0;
           SHK.ρz; SHK.σz; SHK.ρδ; SHK.σδ; SHK.ρs; SHK.σs]
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
    return (; ss = sol.ss, cal, hx = sol.sol_mat.hx, nunstable,
            irf_z = irf(:z, SHK.σz), irf_δ = irf(:δ, SHK.σδ), irf_s = irf(:s, SHK.σs))
end

old = deserialize("R8_snapshot_pre.jls")
new = Dict(k => solve_case(tg) for (k, tg) in cases)
serialize("R8_test_timing.jls", new)

common = intersect(Symbol.(names(old[:p0zero].irf_z)), Symbol.(names(new[:p0zero].irf_z)))
common = setdiff(common, [:δ_e, :u, :N, :v_pret])      # δ_e re-dated; stocks checked below
stocks = [:u, :N, :v_pret]                              # old: states; new: controls

lines = String[]
log!(s) = (push!(lines, s); println(s))
log!("R8 step 2 tests — timing restructure, costless reposting. IRFs in 100 × log deviations, T = $T_IR")

# ── Test 1 ────────────────────────────────────────────────────────────────────
log!("\nTest 1 (p_0 = 0): max |new − old| over t, by shock")
worst_all = 0.0
for shk in (:irf_z, :irf_δ, :irf_s)
    o, n = getfield(old[:p0zero], shk), getfield(new[:p0zero], shk)
    d_common = maximum(maximum(abs.(n[!, v] .- o[!, v])) for v in common)
    d_stocks = maximum(maximum(abs.(n[!, v] .- o[!, v])) for v in stocks)
    d_δe     = maximum(abs.(n.δ_e[2:end] .- o.δ_e[1:end-1]))
    global worst_all = max(worst_all, d_common, d_stocks, d_δe)
    log!(@sprintf("  %-6s common vars %.2e   u,N,v_pret %.2e   δ_e (new[t+1] vs old[t]) %.2e",
                  shk, d_common, d_stocks, d_δe))
end
log!(@sprintf("  Test 1 %s (worst %.2e, tolerance 1e-8)", worst_all < 1e-8 ? "PASSED" : "FAILED", worst_all))

# ── Test 2 ────────────────────────────────────────────────────────────────────
b = new[:base]
log!(@sprintf("\nTest 2 (p_0 = 0.5): solved; eigenvalues of hx with |λ| ≥ 1: %d (want 0)", b.nunstable))
function conservation_gap(r, df)
    lv(v) = getfield(r.ss, v) .* exp.(df[!, v] ./ 100)
    v, L, δe = lv(:v), lv(:L), lv(:δ_e)
    vp, δt  = lv(:v_pret), r.cal.δ .* exp.(df.δ ./ 100)      # δ_t level
    δlag    = r.cal.δ .* exp.(df.δl ./ 100)                  # δ_{t-1} level
    Λ       = (1 .- δe) ./ (1 .- δlag)
    lhs = vp[2:end] .+ L[2:end]
    rhs = Λ[2:end] .* (1 .- δt[1:end-1]) .* (v[1:end-1] .+ L[1:end-1])
    return maximum(abs.(lhs .- rhs))
end
# The identity holds to first order, so in levels its gap is second order in the shock.
# Evaluate it at 1/1000 of a standard deviation (IRFs are linear in the shock size): a
# correct model gives a gap ~1e-6 × the one-SD gap, i.e. at machine-precision scale.
for shk in (:irf_z, :irf_δ, :irf_s)
    df = getfield(b, shk)
    g1, gs = conservation_gap(b, df), conservation_gap(b, df ./ 1000)
    log!(@sprintf("  position identity, %-6s gap in levels: 1 SD %.2e, 1e-3 SD %.2e (ratio %.1e; second order ⇒ ~1e-6)",
                  shk, g1, gs, gs / g1))
end
log!("\n  How much the base IRFs moved vs the pre-R8 code (expected: exit now dated t+1):")
for shk in (:irf_z, :irf_δ, :irf_s), v in (:u, :v, :N, :θ)
    o, n = getfield(old[:base], shk)[!, v], getfield(b, shk)[!, v]
    log!(@sprintf("    %-6s %-3s  old[3,13] = %9.5f %9.5f   new[3,13] = %9.5f %9.5f",
                  shk, v, o[3], o[13], n[3], n[13]))
end
write("R8_test_timing.txt", join(lines, "\n") * "\n")
