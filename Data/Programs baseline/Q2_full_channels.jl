# Q2_full_channels.jl
# =============================================================================
# Exploratory, not draft-bound. TRUE-SHOCK impulse responses (decisions.md S16): which channels
# does the full model's business-formation block move after δ and z shocks, and by how much?
#
# Motivation (Nov 10, 2026): in the data, a one-SD δ shock raises u by about 1.5–2 % with a delayed
# peak at three to four quarters, through depressed job finding (q3_delta_robustness.txt). AGS
# peaks at about 0.4 in the first quarter. Can business formation produce a delayed, larger fall
# in job finding? The variety channel works through w_int = ρ(N) z/μ, ρ = N^(1/(ε−1)), so it needs
# N to move; the discount-factor and user-cost channels work through C and Q.
#
# Cases (as in Q2_persistence_comparison.jl): no exit (p_0 = 0) and the working point (p_0 = 0.8),
# each with priced entry and clone replacement; σ = 1, AGS's reactivation (Λ̄_r 0.25, p_r 0.8,
# ψ_r 1), Block B′'s shock processes. Responses to one-SD shocks, % (100 × log deviation),
# quarterly averages (quarter 0 = the shock's quarter, impact in its first month).
# Variables: N (product lines in production), N_e (entry), ρ (variety), w_int (revenue per worker),
# θ, f (job finding, (1 − η) θ̂ plus the survival term), u, v, e (creation), Q, K, C, and the
# one-quarter-ahead discount factor m (−σ ΔC, quarterly).
#
# Output: Q2_full_channels.txt, Q2_full_channels.png
# =============================================================================

include("run_solution_core.jl")
using Printf, Statistics, CSV, Plots

const SHOCK_CSV = joinpath(@__DIR__, "..", "Bartek analysis", "data", "results", "shock_calibration.csv")
const SHK_B     = shock_processes(SHOCK_CSV; prefix = "bprime")
const REPOST    = (Λr_target = 0.25, p_r = 0.8, ψ_r = 1.0)
const CASES = (
    ("no exit",    merge(TARGETS, (p_0 = 0.0, dest_end_frac = 0.0, Xc_Y = 0.0), REPOST)),
    ("working pt", merge(TARGETS, (p_0 = 0.8,), REPOST)),
)
const VARS = (:N, :N_e, :ρ, :w_int, :θ, :f, :u, :v, :e, :Q, :K, :C, :m)
const NQ = 13                                   # quarters 0–12

lines = String[]
log!(str) = (push!(lines, str); println(str))

function build_model(; clone_N = nothing)
    f = gen_model_equations()
    if !isnothing(clone_N)
        f[4] = Ntp - clone_N
    end
    m = (parameters = parameters, estimate = estimate, estimation = position,
         npar = length(parameters), ns = length(estimate), priors = priors,
         x = x, y = y, xp = xp, yp = yp, variables = variables, varnames = varnames,
         nx = nx, ny = ny, nvar = nvar, e = ex, eta = eta, ne = ne,
         f = f, nf = nvar, SS = SS_symbolics(parameters, TARGETS),
         PAR_SS = parameters[:], flag_order = flag_order, flag_deviation = flag_deviation,
         flag_SSsolver = flag_SSsolver)
    process_model(m)
    return m
end

function pre_exit_N(tg)
    cal = calibrate_shares(tg)
    st = steady_state(cal)
    return st.N / ((1 - st.δ_e) / (1 - cal.δ))
end

quarterly(m) = [mean(m[3k-2:3k]) for k in 1:NQ]

"""Quarterly responses (%, 100 × log deviation) of VARS to a one-SD shock to `state`."""
function responses(mdl, sol, cal, state, sd)
    ir = DataFrame(simulate_model((; mdl..., ne = 1), sol.sol_mat, 3 * NQ + 1, shock_col(state, sd),
                                  sol.SS, true, true), varnames)[2:end, :]
    f_m = log.(1 .- cal.δ .* exp.(ir.δ)) .- log(1 - cal.δ) .+ (1 - cal.η_L) .* ir.θ   # job finding
    out = Dict{Symbol,Vector{Float64}}()
    for v in VARS
        if v == :f
            out[v] = 100 .* quarterly(f_m)
        elseif v == :m
            c = quarterly(ir.C)
            out[v] = 100 .* (-cal.σ) .* vcat(diff(c), 0.0)     # m_{t+1}: −σ (C_{t+1} − C_t), quarterly
        else
            out[v] = 100 .* quarterly(ir[!, v])
        end
    end
    return out
end

function run_case(label, tg, mdl)
    cal = calibrate_shares(tg)
    PAR = build_PAR(cal, SHK_B)
    sol = solution_interface(mdl, PAR, SS_numeric(PAR, tg))
    return Dict(:δ => responses(mdl, sol, cal, :δ, SHK_B.σδ), :z => responses(mdl, sol, cal, :z, SHK_B.σz))
end

res = Dict{String,Any}()
model = build_model()
for (lab, tg) in CASES
    res["$lab, priced"] = run_case("$lab, priced", tg, model)
end
for (lab, tg) in CASES
    res["$lab, clone"] = run_case("$lab, clone", tg, build_model(clone_N = pre_exit_N(tg)))
end
order = ["no exit, priced", "no exit, clone", "working pt, priced", "working pt, clone"]

log!("Full model: channels after δ and z shocks (TRUE shocks; decisions.md S16). σ = 1, AGS reactivation,")
log!(@sprintf("Block B′ shocks: ρz %.4f σz %.4f ρδ %.4f σδ %.4f ρs %.4f σs %.4f (monthly).", SHK_B...))
log!("Quarterly averages of the response to a one-SD shock, % (100 × log deviation); q0 = the shock's quarter.")
log!("Each cell: peak (quarter of the peak) | cumulative over quarters 0–12.")
log!("Data, δ shock (q3_delta_robustness.txt, HP): u peak +1.4 to +1.9 at quarter 3–4; job finding −1.2 to")
log!("−1.6 at quarter 4. AGS (q3_model_svar.txt, through the SVAR): u +0.4 at quarter 1.")
for shk in (:δ, :z)
    log!("\n$(shk) shock")
    log!("  " * rpad("variable", 9) * join([rpad(c, 24) for c in order]))
    for v in VARS
        cells = String[]
        for c in order
            p = res[c][shk][v]
            k = argmax(abs.(p))
            push!(cells, rpad(@sprintf("%+.3f (q%d) | %+.2f", p[k], k - 1, sum(p)), 24))
        end
        log!("  " * rpad(string(v), 9) * join(cells))
    end
end
log!("\nf: gross job finding (1 − δ_t) A θ_t^(1−η). m: −σ times the quarterly change in C (the discount")
log!("factor's log deviation for the next quarter). N: lines in production after exit; ρ = N^(1/(ε−1)).")

open("Q2_full_channels.txt", "w") do io
    foreach(l -> println(io, l), lines)
end

# figure: δ shock (top) and z shock (bottom); N, w_int, f, u. run_solution_core.jl selects the
# PGFPlotsX backend (impulse_response_plots.jl), which needs LaTeX; without it (e.g. the cloud
# container) fall back to GR.
if isnothing(Sys.which("pdflatex"))
    gr()
end
PAL = Dict("no exit, priced" => "#2a78d6", "no exit, clone" => "#eb6834",
           "working pt, priced" => "#1baf7a", "working pt, clone" => "#eda100")
show_vars = (:N, :w_int, :f, :u)
fig = plot(layout = (2, length(show_vars)), size = (1300, 600), background_color = "#fcfcfb",
           titlefontsize = 10, tickfontsize = 7, legendfontsize = 7, gridalpha = 0.15,
           left_margin = 6Plots.mm, bottom_margin = 5Plots.mm)
for (row, shk) in enumerate((:δ, :z)), (col, v) in enumerate(show_vars)
    for c in order
        if row == 1 && col == 1
            lab = c
        else
            lab = ""
        end
        plot!(fig[row, col], 0:NQ-1, res[c][shk][v], lw = 2, color = PAL[c], label = lab)
    end
    hline!(fig[row, col], [0], color = "#8a8984", lw = 1, ls = :dash, label = "")
    title!(fig[row, col], "$(shk) shock: $(v)")
    if row == 2
        xlabel!(fig[row, col], "quarters")
    end
    if col == 1
        ylabel!(fig[row, col], "%")
    end
end
savefig(fig, "Q2_full_channels.png")

