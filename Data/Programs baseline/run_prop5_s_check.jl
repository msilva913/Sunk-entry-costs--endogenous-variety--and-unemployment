# run_prop5_s_check.jl
# =============================================================================
# NUMERICAL CHECK of Proposition 5 Part 1 (prop:ds_asymmetry, app:proof_ds).
#
# Timing convention: the aggregate state, including s_t, is known at the start
# of t, so tightness and entry at t respond to the shock. With p_0 = 0, u_t and
# v_pre,t are predetermined, d = (1-δ_t)(1-u_t), and M_t = f_t u_t:
#   h = 0:  Δu_t = 0,  Δv_t = Δe_t
#   h = 1:  Δu_{t+1}     = d − (1−δ_t)·ΔM_t,        ΔM_t = (1−η_L)·q_t·Δe_t
#           Δv_pre,t+1   = d + (1−δ_t)·(Δe_t − ΔM_t)
#           Δv_{t+1}     = Δv_pre,t+1 + Δe_{t+1}
# Part 1: u and v comove positively at h = 1 provided
#   (i)  (1−δ_t)·ΔM_t / d < 1                         (unemployment rises)
#   (ii) [−(1−δ_t)(Δe_t − ΔM_t) − Δe_{t+1}] / d < 1   (vacancies rise)
# Both ratios are reported; the draft cites their maxima.
#
# DESIGN
#   Benchmark: p_0 = 0, dest_end_frac = 0, Xc_Y = 0 (the exogenous-exit arm of
#   Comparison B), b_ratio = 0.9, x_v = 0.5. σ = 1, DS-CES: the discount-factor
#   and variety channels are ON, so this is the general case.
#   Grid: dest_ann ∈ {0.0320 (D1, BED), 0.0754 (current code)}
#         ρ_s     ∈ {0, 0.5, 0.8741 (part6b), 0.95}
#         ξ_inv   ∈ {0.5, 1, 2} at the part6b ρ_s
#   Shock: one s innovation of size σ_s = 0.0854 (part6b). First order, so the
#   ratios are scale-free.
#
# TIMING (simulate_model): row 1 is the steady state; the s state is shocked in
# row 2 (period t), where θ_t and e_t respond; u_{t+1}, v_pre,t+1, e_{t+1},
# v_{t+1} are in row 3.
# Output: prop5_s_check.csv
# =============================================================================

include("run_solution_core.jl")
using Printf

ρ_z = 0.902
σ_z = 0.0092
ρ_δ = 0.592
σ_δ = 0.0669
σ_s = 0.0854      # part6b
ρ_s_data = 0.8741 # part6b: AR(1) on s = (τ−δ_e)/(1−δ_e), 1992Q3–2019Q4, monthly

T_IR        = 12
flag_IR     = true
flag_logdev = true

bench = (TARGETS...,
    dest_end_frac = 0.0,
    p_0           = 0.0,
    Xc_Y          = 0.0,
    b_ratio       = 0.9,
    x_v           = 0.5)

function build_model_once(dummy_targets)
    f = gen_model_equations()
    SS_sym = SS_symbolics(parameters, dummy_targets)
    model = (parameters = parameters, estimate = estimate, estimation = position,
             npar = length(parameters), ns = length(estimate), priors = priors,
             x = x, y = y, xp = xp, yp = yp, variables = variables, varnames = varnames,
             nx = nx, ny = ny, nvar = nvar, e = ex, eta = eta, ne = ne,
             f = f, nf = nvar, SS = SS_sym, PAR_SS = parameters[:],
             flag_order = flag_order, flag_deviation = flag_deviation,
             flag_SSsolver = flag_SSsolver)
    process_model(model)
    return model
end

function s_check(model, targets_variant, ρ_s_val)
    cal = calibrate_shares(targets_variant)
    @unpack f_e, δ, s, z, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, x_m, ψ, f_m, p_0 = cal
    s < -1e-8 && error("negative s = $s")
    s = max(s, 0.0)
    zbar = z; δbar = δ; sbar = s

    PAR = [f_e; zbar; δbar; sbar; b; ϕ; r; σ; ε; A; η_L; κ; ξ_inv; x_m; ψ; f_m; p_0;
           ρ_z; σ_z; ρ_δ; σ_δ; ρ_s_val; σ_s]
    sol = solution_interface(model, PAR, SS_numeric(PAR, targets_variant))
    @unpack ss, SS, sol_mat = sol

    eta_s_col = reshape([0.0, 0.0, 0.0, 0.0, 0.0, σ_s], nx, 1)
    irf = DataFrame(simulate_model((; model..., ne = 1), sol_mat, T_IR, eta_s_col, SS,
                                   flag_IR, flag_logdev), varnames)

    t, h1 = 2, 3
    surv     = 1 - ss.δ_e
    v_pre_ss = ss.v - ss.e
    d        = surv * (1 - ss.u) * sbar * irf.s[t]

    Δe_t   = ss.e * irf.e[t]
    Δθ_t   = ss.θ * irf.θ[t]
    ΔM_t   = ss.u * (1 - η_L) * ss.q * Δθ_t        # u_t fixed: ΔM = ū·f'(θ̄)·Δθ, f' = (1−η_L)q
    Δu     = ss.u * irf.u[h1]
    Δv_pre = v_pre_ss * irf.v_pret[h1]
    Δe_t1  = ss.e * irf.e[h1]
    Δv     = ss.v * irf.v[h1]

    return (dest_ann = targets_variant.dest_ann, ρ_s = ρ_s_val, ξ_inv = ξ_inv,
            θ_ss = ss.θ, d = d, Δe_t = Δe_t, ΔM_t = ΔM_t, Δe_t1 = Δe_t1,
            Δu = Δu, Δv_pre = Δv_pre, Δv = Δv,
            resid_u    = Δu - (d - surv * ΔM_t),
            resid_vpre = Δv_pre - (d + surv * (Δe_t - ΔM_t)),
            resid_v    = Δv - (Δv_pre + Δe_t1),
            cond_i  = surv * ΔM_t / d,
            cond_ii = (-surv * (Δe_t - ΔM_t) - Δe_t1) / d,
            holds = (Δu > 0) && (Δv > 0))
end

grid = Tuple{Float64,Float64,Float64}[]
for da in (0.0320, 0.0754)
    for ρ in (0.0, 0.5, ρ_s_data, 0.95)
        push!(grid, (da, ρ, 1.0))
    end
    for ξi in (0.5, 2.0)
        push!(grid, (da, ρ_s_data, ξi))
    end
end

println("\nBuilding symbolic model (once) ...")
model = build_model_once((bench..., dest_ann = grid[1][1]))

rows = NamedTuple[]
for (da, ρ, ξi) in grid
    lab = "dest_ann=$da  ρ_s=$ρ  ξ_inv=$ξi"
    try
        push!(rows, s_check(model, (bench..., dest_ann = da, ξ_inv = ξi), ρ))
        println("  ok   ", lab)
    catch err
        println("  FAIL ", lab, "  → ", sprint(showerror, err)[1:min(end, 160)])
    end
end

res = DataFrame(rows)
open("prop5_s_check.csv", "w") do io
    println(io, join(names(res), ","))
    for r in eachrow(res)
        println(io, join(string.(collect(r)), ","))
    end
end

println("\n" * "="^104)
println("Prop. 5 Part 1 (p_0 = 0). Levels ×1e4. Conditions hold iff cond_i < 1 and cond_ii < 1")
println("="^104)
@printf("%-9s %-7s %-6s %-6s %9s %9s %9s %9s %9s %9s %8s %8s %6s\n",
        "dest_ann", "ρ_s", "ξ_inv", "θ̄", "d", "Δe_t", "ΔM_t", "Δe_t+1", "Δu", "Δv",
        "cond_i", "cond_ii", "holds")
for r in eachrow(res)
    @printf("%-9.4f %-7.4f %-6.2f %-6.3f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %8.4f %8.4f %6s\n",
            r.dest_ann, r.ρ_s, r.ξ_inv, r.θ_ss, 1e4*r.d, 1e4*r.Δe_t, 1e4*r.ΔM_t,
            1e4*r.Δe_t1, 1e4*r.Δu, 1e4*r.Δv, r.cond_i, r.cond_ii, string(r.holds))
end
@printf("\nmax cond_i = %.4f   max cond_ii = %.4f   (as shares of the reposting inflow d)\n",
        maximum(res.cond_i), maximum(res.cond_ii))
@printf("identity residuals (should be ≈ 0, relative to d): u %.1e   v_pre %.1e   v %.1e\n",
        maximum(abs.(res.resid_u ./ res.d)), maximum(abs.(res.resid_vpre ./ res.d)),
        maximum(abs.(res.resid_v ./ res.d)))
println("Saved: prop5_s_check.csv")
