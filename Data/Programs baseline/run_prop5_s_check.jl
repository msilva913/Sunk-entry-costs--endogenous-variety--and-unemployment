# run_prop5_s_check.jl
# =============================================================================
# NUMERICAL CHECK of Proposition 5 Part 1 (prop:ds_asymmetry, app:proof_ds).
#
# Part 1 (p_0 = 0: exogenous exit, costless reposting) states that a positive
# s_t shock raises u_{t+1} and pre-committed vacancies v_pre,t+1 by the same
# amount d = (1-δ_t)(1-u_t), so total vacancies rise at h = 1 provided entry
# does not contract by more than the inflow:  -∂e_{t+1}/∂s_t < d.
# The proof signs the entry response analytically only under ρ_s = 0 plus the
# prop:independence conditions (σ = 0, ζ = 0), and then only for free entry or
# θ̄ = 1. This program checks the condition at the calibrated p_0 = 0 benchmark.
#
# DESIGN
#   Benchmark: p_0 = 0, dest_end_frac = 0, Xc_Y = 0 (the exogenous-exit arm of
#   Comparison B), b_ratio = 0.9, x_v = 0.5. σ = 1, DS-CES: the GE channels
#   (discount factor, variety) are ON, so this is the general case the proof
#   leaves to quantitative magnitudes.
#   Grid: dest_ann ∈ {0.0320 (D1, BED), 0.0754 (current code)}
#         ρ_s     ∈ {0, 0.5, 0.8741 (part6b), 0.95}
#         ξ_inv   ∈ {0.5, 1, 2} at the part6b ρ_s
#   Shock: one s innovation of size σ_s = 0.0854 (part6b). First order, so the
#   SIGN of each response does not depend on σ_s; ratios are scale-free.
#
# TIMING (simulate_model): row 1 is the steady state, the s state is shocked in
# row 2 (period t), and u_{t+1}, v_pre,t+1, e_{t+1}, v_{t+1} respond in row 3.
#
# REPORTED (levels, first order: Δx = x̄·x̂)
#   d_theory   (1-δ_e)·(1-ū)·s̄·ŝ_t             analytical inflow (Steps 1-2)
#   Δu, Δv_pre                                  should both equal d_theory
#   Δe, Δv                                      Δv = Δv_pre + Δe
#   ratio  -Δe/Δv_pre                           condition holds iff < 1
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

    t, h1 = 2, 3                       # shock period, response period
    v_pre_ss = ss.v - ss.e
    d_theory = (1 - ss.δ_e) * (1 - ss.u) * sbar * irf.s[t]
    Δu     = ss.u * irf.u[h1]
    Δv_pre = v_pre_ss * irf.v_pret[h1]
    Δe     = ss.e * irf.e[h1]
    Δv     = ss.v * irf.v[h1]
    return (dest_ann = targets_variant.dest_ann, ρ_s = ρ_s_val, ξ_inv = ξ_inv,
            θ_ss = ss.θ, u_ss = ss.u, δ_e_ss = ss.δ_e,
            d_theory = d_theory, Δu = Δu, Δv_pre = Δv_pre, Δe = Δe, Δv = Δv,
            Δv_check = Δv - (Δv_pre + Δe), δ_e_resp = irf.δ_e[h1],
            ratio = -Δe / Δv_pre, holds = Δv > 0)
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

println("\n" * "="^100)
println("Prop. 5 Part 1 check (p_0 = 0 benchmark). Levels ×1e4. Condition: -Δe/Δv_pre < 1 ⇔ Δv > 0")
println("="^100)
@printf("%-9s %-7s %-6s %-6s %9s %9s %9s %9s %9s %8s %6s\n",
        "dest_ann", "ρ_s", "ξ_inv", "θ̄", "d_theory", "Δu", "Δv_pre", "Δe", "Δv", "-Δe/Δvp", "holds")
for r in eachrow(res)
    @printf("%-9.4f %-7.4f %-6.2f %-6.3f %9.4f %9.4f %9.4f %9.4f %9.4f %8.3f %6s\n",
            r.dest_ann, r.ρ_s, r.ξ_inv, r.θ_ss, 1e4*r.d_theory, 1e4*r.Δu, 1e4*r.Δv_pre,
            1e4*r.Δe, 1e4*r.Δv, r.ratio, string(r.holds))
end
println("\nmax |Δv - (Δv_pre + Δe)| = ", maximum(abs.(res.Δv_check)),
        "   max |δ_e response| = ", maximum(abs.(res.δ_e_resp)), "   (both should be ≈ 0)")
println("Saved: prop5_s_check.csv")
