# run_prop5_parts23_check.jl
# =============================================================================
# NUMERICAL CHECKS for Proposition 5 Parts 2 and 3 (prop:ds_asymmetry) under the
# paper's timing convention: the aggregate state (incl. s_t, δ_t) is known at
# the start of t, so θ_t and e_t respond to period-t shocks.
#
# PART 3 (δ shock), p_0 = 0 benchmark. Exact here: with p_0 = 0 the code's exit
# dating issue (model_equations.md "Code status") is irrelevant. u_t, v_pre,t
# predetermined; M_t = f_t u_t; A_u = f̄ū + (1−s̄)(1−ū); A_v = (1−q̄)v̄ + s̄(1−ū):
#   Δu_{t+1}     = A_u·Δδ − (1−δ̄)·ΔM_t
#   Δv_pre,t+1   = −A_v·Δδ + (1−δ̄)·(Δe_t − ΔM_t)
#   Δv_{t+1}     = Δv_pre,t+1 + Δe_{t+1}
# Negative comovement at h = 1 iff
#   c3_u = (1−δ̄)ΔM_t / (A_u Δδ) < 1                       (unemployment rises)
#   c3_v = [(1−δ̄)(Δe_t − ΔM_t) + Δe_{t+1}] / (A_v Δδ) < 1  (vacancies fall)
#
# PART 2 feedback size, p_0 > 0 baseline (PATH B, dest_elast_target = 5). The
# draft's 2×2 Jacobian for (χ^c_{t+1}, N_{t+1}) omits that post-exit
# unemployment u_{t+1} = 1 − F(χ^c_{t+1})·Ẽ_{t+1} depends on the cutoff. With it
# the system is 3×3 and
#   det J3 = 1 − g_N·(1−δ̄)F'(χ̄^c)·B̄ + g_u·F'(χ̄^c)·Ẽ,   g_u < 0,
# so det J3 > 0 needs |g_u|·F'·Ẽ < 1 + |g_N|(1−δ̄)F'B̄ when g_N < 0.
# Also reported: the sign of the period-t exit response Λ̂_t to s and δ shocks.
# ⚠️ INDICATIVE ONLY for p_0 > 0: the code applies δ_e,t (built from Λ_t) to the
# t→t+1 flows, while the draft applies (1−δ_t)Λ_{t+1}. Redo after R8.
# Output: prop5_parts23_check.csv (Part 3 grid)
# =============================================================================

include("run_solution_core.jl")
using Printf

ρ_z = 0.902
σ_z = 0.0092
ρ_δ_data = 0.592
σ_δ = 0.0669
ρ_s = 0.8741
σ_s = 0.0854

T_IR        = 40
flag_IR     = true
flag_logdev = true

bench = (TARGETS..., dest_end_frac = 0.0, p_0 = 0.0, Xc_Y = 0.0, b_ratio = 0.9, x_v = 0.5)
base  = (TARGETS..., dest_end_frac = 0.5, p_0 = 0.5, dest_elast_target = 5.0,
         b_ratio = 0.9, x_v = 0.5, ξ_inv = 1.0)

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

function solve(model, tv, ρ_δ)
    cal = calibrate_shares(tv)
    @unpack f_e, δ, s, z, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, x_m, ψ, f_m, p_0 = cal
    s = max(s, 0.0)
    PAR = [f_e; z; δ; s; b; ϕ; r; σ; ε; A; η_L; κ; ξ_inv; x_m; ψ; f_m; p_0;
           ρ_z; σ_z; ρ_δ; σ_δ; ρ_s; σ_s]
    sol = solution_interface(model, PAR, SS_numeric(PAR, tv))
    return cal, sol
end

function irf_of(model, sol, col)
    @unpack SS, sol_mat = sol
    eta_col = reshape(col, nx, 1)
    DataFrame(simulate_model((; model..., ne = 1), sol_mat, T_IR, eta_col, SS,
                             flag_IR, flag_logdev), varnames)
end

function part3_check(model, tv, ρ_δ)
    cal, sol = solve(model, tv, ρ_δ)
    ss = sol.ss
    irf = irf_of(model, sol, [0.0, 0.0, 0.0, 0.0, σ_δ, 0.0])
    t, h1 = 2, 3
    δbar, sbar, η_L = cal.δ, cal.s, cal.η_L
    surv = 1 - ss.δ_e
    f̄ = ss.θ * ss.q
    Δδ = δbar * irf.δ[t]
    A_u = f̄ * ss.u + (1 - sbar) * (1 - ss.u)
    A_v = (1 - ss.q) * ss.v + sbar * (1 - ss.u)
    Δe_t   = ss.e * irf.e[t]
    ΔM_t   = ss.u * (1 - η_L) * ss.q * ss.θ * irf.θ[t]
    Δu     = ss.u * irf.u[h1]
    Δv_pre = (ss.v - ss.e) * irf.v_pret[h1]
    Δe_t1  = ss.e * irf.e[h1]
    Δv     = ss.v * irf.v[h1]
    return (dest_ann = tv.dest_ann, ρ_δ = ρ_δ, ξ_inv = cal.ξ_inv,
            Δu = Δu, Δv = Δv, Δe_t = Δe_t, Δe_t1 = Δe_t1, ΔM_t = ΔM_t,
            resid_u = (Δu - (A_u * Δδ - surv * ΔM_t)) / (A_u * Δδ),
            resid_vpre = (Δv_pre - (-A_v * Δδ + surv * (Δe_t - ΔM_t))) / (A_v * Δδ),
            c3_u = surv * ΔM_t / (A_u * Δδ),
            c3_v = (surv * (Δe_t - ΔM_t) + Δe_t1) / (A_v * Δδ),
            holds = (Δu > 0) && (Δv < 0))
end

grid = Tuple{Float64,Float64,Float64}[]
for da in (0.0320, 0.0754)
    for ρ in (0.0, ρ_δ_data, 0.9)
        push!(grid, (da, ρ, 1.0))
    end
    for ξi in (0.5, 2.0)
        push!(grid, (da, ρ_δ_data, ξi))
    end
end

println("\nBuilding symbolic model (once) ...")
model = build_model_once((bench..., dest_ann = 0.0320))

rows = NamedTuple[]
for (da, ρ, ξi) in grid
    try
        push!(rows, part3_check(model, (bench..., dest_ann = da, ξ_inv = ξi), ρ))
    catch err
        println("  FAIL dest_ann=$da ρ_δ=$ρ ξ_inv=$ξi → ", sprint(showerror, err)[1:min(end, 160)])
    end
end
res = DataFrame(rows)
open("prop5_parts23_check.csv", "w") do io
    println(io, join(names(res), ","))
    for r in eachrow(res)
        println(io, join(string.(collect(r)), ","))
    end
end

println("\n" * "="^96)
println("PART 3 (δ shock, p_0 = 0). Levels ×1e4. Negative comovement iff c3_u < 1 and c3_v < 1")
println("="^96)
@printf("%-9s %-6s %-6s %9s %9s %9s %9s %9s %8s %8s %6s\n",
        "dest_ann", "ρ_δ", "ξ_inv", "Δu", "Δv", "Δe_t", "Δe_t+1", "ΔM_t", "c3_u", "c3_v", "holds")
for r in eachrow(res)
    @printf("%-9.4f %-6.3f %-6.2f %9.4f %9.4f %9.4f %9.4f %9.4f %8.4f %8.4f %6s\n",
            r.dest_ann, r.ρ_δ, r.ξ_inv, 1e4*r.Δu, 1e4*r.Δv, 1e4*r.Δe_t, 1e4*r.Δe_t1,
            1e4*r.ΔM_t, r.c3_u, r.c3_v, string(r.holds))
end
@printf("max c3_u = %.4f   max c3_v = %.4f   max |identity residual| u %.1e  v_pre %.1e\n",
        maximum(res.c3_u), maximum(res.c3_v),
        maximum(abs.(res.resid_u)), maximum(abs.(res.resid_vpre)))

# ── Part 3 beyond h = 1: vacancy and unemployment paths after a δ shock (p_0 = 0) ──
println("\nPART 3 paths (p_0 = 0, ρ_δ = $ρ_δ_data, ξ_inv = 1). Levels ×1e4 by month after the shock")
for da in (0.0320, 0.0754)
    cal, sol = solve(model, (bench..., dest_ann = da), ρ_δ_data)
    ss = sol.ss
    irf = irf_of(model, sol, [0.0, 0.0, 0.0, 0.0, σ_δ, 0.0])
    hs = 2:min(nrow(irf), 26)
    println("  dest_ann=$da  Δv: ", join([@sprintf("%+.3f", 1e4 * ss.v * irf.v[h]) for h in hs], " "))
    println("  dest_ann=$da  Δu: ", join([@sprintf("%+.3f", 1e4 * ss.u * irf.u[h]) for h in hs], " "))
end

# ── Part 3 at quarterly frequency: every spec in the grid ─────────────────────
# Quarter q averages months 3q..3q+2 counted from the shock month (row 2), the
# frequency of the LPs and of estimation. Reports the largest quarterly-average
# Δv and the smallest quarterly-average Δu over the first 12 quarters, and the
# first month from which Δv stays negative.
qrows = NamedTuple[]
for (da, ρ, ξi) in grid
    cal, sol = solve(model, (bench..., dest_ann = da, ξ_inv = ξi), ρ)
    ss = sol.ss
    irf = irf_of(model, sol, [0.0, 0.0, 0.0, 0.0, σ_δ, 0.0])
    v = [ss.v * irf.v[h] for h in 2:nrow(irf)]
    u = [ss.u * irf.u[h] for h in 2:nrow(irf)]
    nq = div(length(v), 3)
    vq = [sum(v[3q+1:3q+3]) / 3 for q in 0:nq-1]
    uq = [sum(u[3q+1:3q+3]) / 3 for q in 0:nq-1]
    first_neg = findfirst(m -> all(v[m:end] .< 0), 1:length(v))
    push!(qrows, (dest_ann = da, ρ_δ = ρ, ξ_inv = ξi, v_h1 = v[2],
                  max_vq = maximum(vq), min_uq = minimum(uq),
                  first_neg_month = isnothing(first_neg) ? -1 : first_neg - 1))
end
qres = DataFrame(qrows)
open("prop5_part3_quarterly.csv", "w") do io
    println(io, join(names(qres), ","))
    for r in eachrow(qres)
        println(io, join(string.(collect(r)), ","))
    end
end
println("\nPART 3 at quarterly frequency (p_0 = 0). Levels ×1e4; quarters average months 3q..3q+2")
@printf("%-9s %-6s %-6s %9s %10s %10s %16s\n", "dest_ann", "ρ_δ", "ξ_inv", "Δv(h=1)",
        "max Δv_q", "min Δu_q", "Δv<0 from month")
for r in eachrow(qres)
    @printf("%-9.4f %-6.3f %-6.2f %9.4f %10.4f %10.4f %16d\n", r.dest_ann, r.ρ_δ, r.ξ_inv,
            1e4 * r.v_h1, 1e4 * r.max_vq, 1e4 * r.min_uq, r.first_neg_month)
end

# ── Part 2 feedback size and period-t exit sign, p_0 > 0 baseline (indicative) ──
println("\n" * "="^96)
println("PART 2 diagnostics at the p_0 = 0.5 baseline (PATH B). INDICATIVE: code exit dating ≠ draft")
println("="^96)
for da in (0.0320, 0.0754)
    try
        cal, sol = solve(model, (base..., dest_ann = da), ρ_δ_data)
        ss = sol.ss
        @unpack p_0, ψ, f_m, ε, f_e = cal
        zl = cal.z
        μ = ε / (ε - 1)
        Λ = (1 - ss.δ_e) / (1 - cal.δ)
        Fp = ss.x_c < f_m ? p_0 * ψ * (ss.x_c / f_m)^ψ / ss.x_c : 0.0
        g(N, u) = N^(1/(ε-1)) * zl * (1 - u) * (μ - 1) / (μ * N) + f_e * N^(1/(ε-1)) / μ
        hN, hu = 1e-6 * ss.N, 1e-6
        g_N = (g(ss.N + hN, ss.u) - g(ss.N - hN, ss.u)) / (2hN)
        g_u = (g(ss.N, ss.u + hu) - g(ss.N, ss.u - hu)) / (2hu)
        B = ss.N + ss.N_e
        Ẽ = (1 - ss.u) / Λ
        termN = -g_N * (1 - cal.δ) * Fp * B
        termu = g_u * Fp * Ẽ
        irf_s = irf_of(model, sol, [0.0, 0.0, 0.0, 0.0, 0.0, σ_s])
        irf_d = irf_of(model, sol, [0.0, 0.0, 0.0, 0.0, σ_δ, 0.0])
        Λhat(irf, row) = log(1 - ss.δ_e * exp(irf.δ_e[row] / 1)) - log(1 - ss.δ_e) -
                         (log(1 - cal.δ * exp(irf.δ[row])) - log(1 - cal.δ))
        @printf("dest_ann=%.4f  Λ̄=%.5f  F'(χ^c)=%.3e  g_N=%.3e  g_u=%.3e\n", da, Λ, Fp, g_N, g_u)
        @printf("   det J3 = 1 %+.2e (N loop) %+.2e (u loop) = %.6f   [old 2x2 det = %.6f]\n",
                termN, termu, 1 + termN + termu, 1 + termN)
        @printf("   period-t exit response Λ̂_t: s shock %+.3e   δ shock %+.3e   (log dev.)\n",
                Λhat(irf_s, 2), Λhat(irf_d, 2))
        @printf("   x̂_c,t:  s shock %+.3e   δ shock %+.3e\n", irf_s.x_c[2], irf_d.x_c[2])
        # Role of endogenous exit (role_of_ingredients.md): cutoff composition and
        # the endogenous-exit response relative to the exogenous δ impulse.
        irf_z = irf_of(model, sol, [0.0, 0.0, 0.0, σ_z, 0.0, 0.0])
        profit = ss.Y_c * (μ - 1) / (μ * ss.N)
        ς = (Λ - (1 - p_0)) / p_0
        el = ψ * ς / (1 - ς)                      # elasticity of the endogenous exit rate to χ^c
        x_end = 1 - Λ                             # endogenous exit rate
        @printf("   cutoff χ^c = %.4f: per-line profit %.4f (%.1f%%), option value ν_f %.4f (%.1f%%)\n",
                ss.x_c, profit, 100 * profit / ss.x_c, ss.ν_f, 100 * ss.ν_f / ss.x_c)
        @printf("   endogenous exit rate %.5f/month vs τ %.4f; exit elasticity to χ^c = %.2f\n",
                x_end, ss.δ_e + cal.s * (1 - ss.δ_e), el)
        for (lab, ir, sz) in (("z", irf_z, σ_z), ("δ", irf_d, σ_δ), ("s", irf_s, σ_s))
            pk = maximum(abs.(ir.x_c[2:end]))
            Δx = el * x_end * pk                  # peak level change in the endogenous exit rate
            @printf("   %s shock: peak |x̂_c| %.2e → peak Δ(endog. exit) %.2e/month; δ impulse δ̄·σ_δ = %.2e\n",
                    lab, pk, Δx, cal.δ * σ_δ)
        end
        # Model counterpart of the part6b VAR coefficient â_21 (exit on lagged z, log
        # cycles, quarterly; data −1.05, HC3 SE 0.83): peak response of log δ_e to a
        # 1% fall in z. Rough: one-quarter VAR coefficient vs. a peak IRF ratio.
        el_z = maximum(abs.(irf_z.δ_e[2:end])) / σ_z   # IRFs here are log deviations (not ×100)
        @printf("   elasticity of δ_e to z (peak |δ̂_e| per 1%% z): %.3f   [data â_21 ≈ -1.05, SE 0.83]\n", el_z)
    catch err
        println("  FAIL baseline dest_ann=$da → ", sprint(showerror, err)[1:min(end, 200)])
    end
end
println("\nSaved: prop5_parts23_check.csv")
