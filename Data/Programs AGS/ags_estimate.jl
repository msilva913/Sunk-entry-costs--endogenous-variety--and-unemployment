# ags_estimate.jl
# =============================================================================
# Indirect-inference estimation of AGS on the z and s impulse responses (decisions.md D18).
#
# Auxiliary system: SVAR (z, s, δ, u, v, f), recursive in that order, one-SD shocks, h = 0–12, on
# data-length windows filtered as the data are (HP or Hamilton; ENV AGS_FILTER = "hp" or "ham").
# Data targets and SDs: Bartek analysis/data/results/estimation_targets.csv (q3_estimation_targets.py).
#   Loss L(θ) = Σ ((m̄(θ) − d)/sd)² over the 78 responses of u, v, f to z and s (diagonal weights,
#   D18 (a)); m̄ is the mean over R windows. The 39 responses to δ are held out (the test).
# Parameters θ: σ ∈ [0, 4], Λ̄_r ∈ [0.03, 0.99], ξ ∈ [0.25, 10] (ξ_inv = 1/ξ), x_v ∈ [0.05, 0.9],
#   b/w ∈ [0.5, 0.97]; η_L = 0.6 fixed; the other long-run targets met exactly by calibrate_ags.
#   Costly reactivation with p_r = 0.99, as ags_reactivation_grid.jl.
# Shock processes: matched per filter on the SVAR window (shock_calibration_window.csv), not
#   estimated. simulate_model reseeds on every call: common random numbers across θ.
#
# Modes (ENV AGS_MODE):
#   check   the Julia SVAR on the data cycles against Python's targets; one loss evaluation.
#   search  coarse grid over (σ, Λ̄_r), then Nelder–Mead from the two best points (R = AGS_R,
#           default 50); every evaluation is printed and appended to the log.
#   final   at the best logged point: R = 200 evaluation, the δ test (J_δ against its simulated
#           distribution), the Jacobian of the moments for the sandwich, and the model responses.
#
# Outputs (Programs AGS): ags_estimate_log_<filter>.csv (every evaluation), and in final mode
#   ags_estimate_final_<filter>.csv (θ̂, J), ags_estimate_irf_<filter>.csv (model and data
#   responses), ags_estimate_jacobian_<filter>.csv, ags_estimate_jdelta_<filter>.csv.
# =============================================================================

include(joinpath(@__DIR__, "ags_model.jl"))
using CSV, Optim, Printf, Statistics, LinearAlgebra

const FILTER = get(ENV, "AGS_FILTER", "hp")
const MODE   = get(ENV, "AGS_MODE", "check")
const R_SEARCH = parse(Int, get(ENV, "AGS_R", "50"))
const RESULTS = joinpath(@__DIR__, "..", "Bartek analysis", "data", "results")
const NAMES = ["z", "s", "delta", "u", "v", "jf"]
const RESP  = ["u", "v", "jf"]
const H_Q = 12
const BURN = 40
const N_FILTER, N_KEEP = 276, 110                         # data: 1951Q1–2019Q4; BED from 1992Q3
const LOG = joinpath(@__DIR__, "ags_estimate_log_$(FILTER).csv")

# ── targets ─────────────────────────────────────────────────────────────────
const TG_DATA = let d = CSV.read(joinpath(RESULTS, "estimation_targets.csv"), DataFrame)
    d[d.filter .== FILTER, :]
end
const P_LAG = TG_DATA.p[1]
"""Moment order: shock (z, s | δ), response (u, v, jf), horizon 0–12 — as q3_estimation_targets.py."""
function target_vec(block)
    rows = TG_DATA[TG_DATA.block .== block, :]
    return rows.data, rows.sd
end
const D_EST, SD_EST = target_vec("estimation")
const D_TST, SD_TST = target_vec("test")

# ── SVAR, as q3_data_svar.py (OLS with constant, Σ = U'U/(T − k), Cholesky in NAMES order) ──
function var_irf(Y::Matrix{Float64}, p::Int, H::Int)
    T, k = size(Y)
    X = hcat(ones(T - p), [Y[p-j+1:T-j, :] for j in 1:p]...)
    Yt = Y[p+1:end, :]
    B = X \ Yt
    U = Yt - X * B
    Σ = U' * U / (size(U, 1) - size(X, 2))
    A = [Matrix(B[2+(j-1)*k:1+j*k, :]') for j in 1:p]
    P = Matrix(cholesky(Hermitian(Σ)).L)
    Psi = [Matrix{Float64}(I, k, k)]
    for h in 1:H
        M = zeros(k, k)
        for j in 1:min(h, p)
            M += A[j] * Psi[h-j+1]
        end
        push!(Psi, M)
    end
    return [Psi[h+1] * P for h in 0:H]                  # Θ_h, (var, shock)
end

function moment_vec(Θ, shocks)
    out = Float64[]
    for s in shocks
        j = findfirst(==(s), NAMES)
        for r in RESP
            i = findfirst(==(r), NAMES)
            for h in 0:H_Q
                push!(out, Θ[h+1][i, j])
            end
        end
    end
    return out
end

# ── model side ──────────────────────────────────────────────────────────────
const SHK_F = shock_processes(joinpath(RESULTS, "shock_calibration_window.csv"); prefix = FILTER)
const SDS_F = (SHK_F.σz, SHK_F.σδ, SHK_F.σs)
function filter_fn(name)
    if name == "hp"
        return x -> hp_filter(x, LAM)
    else
        return x -> hamilton_filter(x)
    end
end
const FILTER_FN = filter_fn(FILTER)
const LO = [0.0, 0.03, 0.25, 0.05, 0.5]
const HI = [4.0, 0.99, 10.0, 0.9, 0.97]
const PNAMES = ["sigma", "lambda_r", "xi", "x_v", "b_ratio"]

to_unit(θ) = [log((θ[i] - LO[i]) / (HI[i] - θ[i])) for i in 1:5]          # bounds → ℝ
from_unit(u) = [LO[i] + (HI[i] - LO[i]) / (1 + exp(-u[i])) for i in 1:5]

targets_of(θ) = merge(AGS_TARGETS, (σ = θ[1], Λr_target = θ[2], p_r = 0.99, ξ_inv = 1 / θ[3],
                                    x_v = θ[4], b_ratio = θ[5]))

"""Per-window responses (vector of Θ lists) for R data-length windows at θ."""
function window_irfs(θ, R)
    r = solve_ags(targets_of(θ); shk = SHK_F)
    eta_m = hcat([shock_col(i, SDS_F[i]) for i in 1:3]...)
    T_m = 3 * R * (N_FILTER + BURN)
    sim = DataFrame(simulate_model((; model..., ne = 3), r.sol_mat, T_m, eta_m, r.SS, false, true),
                    varnames)
    q = model_quarterly(sim, r.cal)
    out = Vector{Vector{Matrix{Float64}}}()
    for k in 1:R
        last_q = k * (N_FILTER + BURN)
        w = window_filter(q[last_q-N_FILTER+1:last_q, :], N_FILTER, N_KEEP; short = ["δ"],
                          filter = FILTER_FN)
        w = rename(w[!, ["z", "s", "δ", "u", "v", "f"]], "δ" => "delta", "f" => "jf")
        Y = Matrix{Float64}(dropmissing(w))
        push!(out, var_irf(Y, P_LAG, H_Q))
    end
    return out, r
end

mean_irf(list) = [mean([Θ[h] for Θ in list]) for h in 1:H_Q+1]

const BEST = Ref(Inf)
function loss(θ; R = R_SEARCH, stage = "")
    t0 = time()
    status, Jz, Js = "ok", NaN, NaN
    try
        list, _ = window_irfs(θ, R)
        m = moment_vec(mean_irf(list), ("z", "s"))
        e2 = ((m .- D_EST) ./ SD_EST) .^ 2
        Jz, Js = sum(e2[1:39]), sum(e2[40:78])
    catch err
        status = first(split(sprint(showerror, err), '\n'))
    end
    if isnan(Jz)
        J = 1e6
    else
        J = Jz + Js
    end
    BEST[] = min(BEST[], J)
    println(@sprintf("[%s] σ %.3f Λ %.3f ξ %.3f x_v %.3f b/w %.3f | J_z %7.1f J_s %7.1f J %8.1f | best %8.1f | %.1fs %s",
                     stage, θ..., Jz, Js, J, BEST[], time() - t0, status))
    row = DataFrame(stage = stage, sigma = θ[1], lambda_r = θ[2], xi = θ[3], x_v = θ[4],
                    b_ratio = θ[5], R = R, J_z = Jz, J_s = Js, J = J, status = status)
    CSV.write(LOG, row; append = isfile(LOG))
    return J
end

# ── modes ───────────────────────────────────────────────────────────────────
"""The Julia SVAR on the HP data cycles against Python's targets (HP only: the Hamilton cycles are
not written to disk), then one loss evaluation at the default calibration."""
function check()
    cyc = CSV.read(joinpath(RESULTS, "q3_data_svar_cycles_hp.csv"), DataFrame; missingstring = "")
    start = TG_DATA.start[1]
    cyc = cyc[findfirst(==(start), string.(cyc.quarter)):end, :]
    Y = Matrix{Float64}(dropmissing(cyc[!, NAMES]))
    Θ = var_irf(Y, P_LAG, H_Q)
    gap = maximum(abs.(moment_vec(Θ, ("z", "s")) .- D_EST))
    println(@sprintf("Check: Julia SVAR on the data, T = %d, p = %d; max |Julia − Python| = %.2e", size(Y, 1), P_LAG, gap))
    loss([1.0, 0.25, 1.0, 0.5, 0.9]; stage = "check")
end

function search()
    pts = Vector{Vector{Float64}}()
    for σ in (0.0001, 1.0, 2.0, 4.0 - 1e-4), Λ in (0.031, 0.25, 0.9)
        θ = [σ, Λ, 1.0, 0.5, 0.9]
        push!(pts, θ)
        loss(θ; stage = "grid")
    end
    log = CSV.read(LOG, DataFrame)
    g = sort(log[(log.stage .== "grid") .& (log.status .== "ok"), :], :J)
    for k in 1:min(2, nrow(g))
        θ0 = [g.sigma[k], g.lambda_r[k], g.xi[k], g.x_v[k], g.b_ratio[k]]
        θ0 = clamp.(θ0, LO .+ 1e-3, HI .- 1e-3)
        res = optimize(u -> loss(from_unit(u); stage = "nm$k"), to_unit(θ0), NelderMead(),
                       Optim.Options(iterations = 150, show_trace = false))
        println("Nelder–Mead $k: ", from_unit(Optim.minimizer(res)), "  J = ", Optim.minimum(res))
    end
end

function final()
    log = CSV.read(LOG, DataFrame)
    ok = log[(log.status .== "ok") .& (log.R .== R_SEARCH), :]
    b = ok[argmin(ok.J), :]
    θ = [b.sigma, b.lambda_r, b.xi, b.x_v, b.b_ratio]
    R = 200
    list, r = window_irfs(θ, R)
    Θm = mean_irf(list)
    m_est, m_tst = moment_vec(Θm, ("z", "s")), moment_vec(Θm, ("delta",))
    e2 = ((m_est .- D_EST) ./ SD_EST) .^ 2
    Jd = sum(((m_tst .- D_TST) ./ SD_TST) .^ 2)
    # δ test: J_δ of each window against the model mean (its distribution if the model is true)
    Jw = [sum(((moment_vec(Θ, ("delta",)) .- m_tst) ./ SD_TST) .^ 2) for Θ in list]
    pval = mean(Jw .>= Jd)
    CSV.write(joinpath(@__DIR__, "ags_estimate_final_$(FILTER).csv"),
              DataFrame(param = vcat(PNAMES, ["phi", "J_z", "J_s", "J_delta", "p_delta", "R"]),
                        value = vcat(θ, [r.cal.ϕ, sum(e2[1:39]), sum(e2[40:78]), Jd, pval, R])))
    CSV.write(joinpath(@__DIR__, "ags_estimate_jdelta_$(FILTER).csv"), DataFrame(J_window = Jw))
    irf_rows = DataFrame()
    labels = [(s, r_, h) for s in ("z", "s", "delta") for r_ in RESP for h in 0:H_Q]
    mv = vcat(m_est, m_tst)
    dv, sv = vcat(D_EST, D_TST), vcat(SD_EST, SD_TST)
    for (n, (s, r_, h)) in enumerate(labels)
        push!(irf_rows, (shock = s, response = r_, h = h, model = mv[n], data = dv[n], sd = sv[n]))
    end
    CSV.write(joinpath(@__DIR__, "ags_estimate_irf_$(FILTER).csv"), irf_rows)
    # Jacobian of the estimation moments (central differences, 1 % of each range)
    G = zeros(78, 5)
    for i in 1:5
        δθ = 0.01 * (HI[i] - LO[i])
        up, dn = copy(θ), copy(θ)
        up[i] = min(θ[i] + δθ, HI[i] - 1e-6)
        dn[i] = max(θ[i] - δθ, LO[i] + 1e-6)
        mu = moment_vec(mean_irf(window_irfs(up, R)[1]), ("z", "s"))
        md = moment_vec(mean_irf(window_irfs(dn, R)[1]), ("z", "s"))
        G[:, i] = (mu .- md) ./ (up[i] - dn[i])
    end
    CSV.write(joinpath(@__DIR__, "ags_estimate_jacobian_$(FILTER).csv"), DataFrame(G, PNAMES))
    println(@sprintf("Final (R = %d): θ = %s; J_z %.1f J_s %.1f; J_δ %.1f, p = %.3f",
                     R, string(round.(θ; digits = 3)), sum(e2[1:39]), sum(e2[40:78]), Jd, pval))
end

if MODE == "check"
    check()
elseif MODE == "search"
    search()
elseif MODE == "final"
    final()
else
    error("AGS_MODE must be check, search or final")
end
