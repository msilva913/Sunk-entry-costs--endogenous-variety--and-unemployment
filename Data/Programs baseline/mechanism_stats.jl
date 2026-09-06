# mechanism_stats.jl
# =============================================================================
# Draft-bound statistics for Section 5.3 ("Inspecting the mechanism") and the
# Comparison D appendix.
#
# WHY THIS FILE EXISTS
# --------------------
# principles.md N15: every number that appears in the manuscript must be
# produced by a program in the replication workflow, never computed in an
# ad-hoc session and pasted in. The mechanism runners serialize full IRF paths
# and the plotters draw them, but neither emits the scalars the prose quotes
# (half-lives, peak ratios, entry costs, the exit-flow decomposition). This
# script closes that gap: it reads the same `irf_*.jls` the figures are drawn
# from and writes every quoted quantity to disk.
#
# Run AFTER the four mechanism runners, since it consumes their output:
#   run_solution_all_delta.jl        -> irf_all_delta.jls    (Comparison A)
#   run_solution_endog_exit.jl       -> irf_endog_exit.jls   (Comparison B)
#   run_solution_no_variety.jl       -> irf_no_variety.jls   (Comparison C)
#   run_solution_entry_elasticity.jl -> irf_xi_inv.jls       (Comparison D)
#
# OUTPUT
#   mechanism_stats.txt  human-readable, grouped by where the number appears
#   mechanism_stats.tex  \newcommand macros, so the draft can \input them and
#                        stop carrying hard-coded figures altogether
#
# SIGN CONVENTION
#   The z-shock figures plot a NEGATIVE technology shock so that unemployment
#   rises in both panels (see plot_*.jl, flip_sign=true). The serialized IRFs
#   are for a POSITIVE shock, so z-shock quantities are negated here to match
#   what the draft describes. Destruction-shock IRFs are used as stored.
# =============================================================================

using Serialization, DataFrames, Printf

cd(@__DIR__)

const TAU_M = 0.031   # monthly total separation rate (TARGETS.sep); the
                      # denominator for every δ_e/τ ratio the draft reports

# ── Helpers ──────────────────────────────────────────────────────────────────

"""
    to_pp(logdev, ss_level)

Convert a 100×log-deviation path to a percentage-point deviation of a rate.
Unemployment and vacancy rates are quoted in pp in the draft (principles.md §2),
while stocks and asset values are quoted in log deviations.
"""
to_pp(logdev, ss_level) = (exp.(logdev ./ 100) .- 1) .* ss_level .* 100

"""
    half_life(rho)

Months for an AR(1) deviation to decay by half. Applied both to the exogenous
shocks and to the firm stock, whose own autoregressive root is (1-δ_e) from the
N law of motion. Comparing the two is the point: the firm stock is the slow
clock that gives the model its persistence, the shocks are the fast one.
"""
half_life(rho) = log(0.5) / log(rho)

"""
    peak(path)

Largest absolute deviation and the month it occurs. The draft quotes peaks
rather than impact responses because the interesting dynamics here are
propagation through the firm stock, which is predetermined and so contributes
nothing on impact.
"""
function peak(path)
    i = argmax(abs.(path))
    return (value = path[i], month = i - 1)
end

"""
    decay_half_life(path)

Months from the peak until the response falls below half its peak. This is the
model's analogue of IRF persistence and is what the draft means when it says one
specification "recovers faster" than another. Returns -1 when the response has
not halved within the simulated horizon, which is itself informative.
"""
function decay_half_life(path)
    i = argmax(abs.(path))
    j = findfirst(k -> abs(path[k]) < 0.5 * abs(path[i]), i:length(path))
    return j === nothing ? -1 : j - 1
end

# Pull a variable from a serialized case, in the units the draft quotes it.
series(case, shock, var; pp = false, flip = false) = begin
    raw = getfield(case, shock)[!, var]
    x   = pp ? to_pp(raw, getfield(case.ss, var)) : raw
    flip ? -x : x
end

# Accumulates every reported number so both output files stay in sync.
const ROWS = Tuple{String,String,Float64,String}[]
record!(section, key, value, note) = push!(ROWS, (section, key, Float64(value), note))

# ── Load ─────────────────────────────────────────────────────────────────────

A = deserialize("irf_all_delta.jls")      # base vs high_δ
B = deserialize("irf_endog_exit.jls")     # endog vs exog
C = deserialize("irf_no_variety.jls")     # base vs novar
D = deserialize("irf_xi_inv.jls")         # base vs low_ξ

# The mechanism-section baseline is the endogenous-exit calibration, shared by
# Comparisons B, C and D. Comparison A deliberately shuts endogenous exit off
# (p_0 = 0), so its "base" is a different economy and is handled separately.
mech = B.endog

# =============================================================================
# 1. Section opening: steady state and the two clocks
# =============================================================================
# The section opens by contrasting how long a firm lives against how long a
# shock lasts. That gap is the reason a shock which dies within a month can move
# unemployment for years: persistence is inherited from the firm stock, not from
# the driving process.

δe   = mech.ss.δ_e
ρ_z, ρ_δ = 0.902, 0.592     # monthly persistences from part6b (see data_and_files.md)

record!("opening", "SSdeltae",      100δe,                   "monthly destruction rate, %")
record!("opening", "SSdeltabar",    100 * 0.5 * δe,          "exogenous component (ω_δ=0.5), %")
record!("opening", "SStau",         100TAU_M,                "total separation rate, %")
record!("opening", "SSu",           100 * mech.ss.u,         "unemployment rate, %")
record!("opening", "FirmLifeMonths", 1 / δe,                 "expected firm life 1/δ_e, months")
record!("opening", "FirmHalfLife",  half_life(1 - δe),       "half-life of an N deviation, months")
record!("opening", "FirmHalfLifeYears", half_life(1 - δe)/12, "same, years")
record!("opening", "HalfLifeDelta", half_life(ρ_δ),          "δ shock half-life, months")
record!("opening", "HalfLifeZ",     half_life(ρ_z),          "z shock half-life, months")
record!("opening", "ClockRatio",    half_life(1 - δe)/half_life(ρ_δ),
        "how many times slower the firm stock is than the δ shock")
record!("opening", "DeltaTauRatio", δe / TAU_M,              "δ_e/τ, the CK-comparison ratio")
record!("opening", "ReplacementRate", 100δe / (1 - δe),      "N^e/N = δ_e/(1-δ_e), %")

# =============================================================================
# 2. Comparison A: the level of destruction
# =============================================================================
# A raises δ_e from the calibrated value to the Coles-Kelishomi convention
# δ_e = τ. Two forces push the same way. The LOM loading on entry equals δ_e, so
# a high-δ_e firm stock rebuilds faster. And a vacancy attached to a
# short-lived product line is a short-duration asset, so its price tracks
# current conditions more closely. The K-vs-Q contrast below isolates the
# second: the dividend moves almost identically across the two, the asset price
# does not.

δe_A_base, δe_A_high = A.base.ss.δ_e, A.high_δ.ss.δ_e
record!("compA", "AdeltaeBase",  100δe_A_base,               "baseline δ_e, %/month")
record!("compA", "AdeltaeHigh",  100δe_A_high,               "high-δ δ_e, %/month")
record!("compA", "AdeltaeRatio", δe_A_high / δe_A_base,      "how many times larger")
record!("compA", "ASurvivalBase", 1 - δe_A_base,             "(1-δ_e), the N survival coefficient")
record!("compA", "AFirmLifeBase", 1 / δe_A_base / 12,        "expected vacancy/firm life, years")
record!("compA", "AFirmLifeHigh", 1 / δe_A_high / 12,        "same, high-δ, years")

for (tag, case) in (("Base", A.base), ("High", A.high_δ))
    # z shock (plotted as negative, hence flip)
    record!("compA", "Az$(tag)K",  peak(series(case, :irf_z, :K;   flip=true)).value, "z: K̂ peak, %")
    record!("compA", "Az$(tag)Q",  peak(series(case, :irf_z, :Q;   flip=true)).value, "z: Q̂ peak, %")
    record!("compA", "Az$(tag)Ne", peak(series(case, :irf_z, :N_e; flip=true)).value, "z: N̂^e peak, %")
    # δ shock
    p_u = peak(series(case, :irf_δ, :u; pp=true))
    record!("compA", "Ad$(tag)N",   peak(series(case, :irf_δ, :N)).value,  "δ: N̂ trough, %")
    record!("compA", "Ad$(tag)U",   p_u.value,                             "δ: u peak, pp")
    record!("compA", "Ad$(tag)Uhl", decay_half_life(series(case, :irf_δ, :u; pp=true)),
            "δ: months for u to fall to half its peak")
end

# =============================================================================
# 3. Comparison B: endogenous versus exogenous exit
# =============================================================================
# B removes the continuation-cost margin while holding δ̄ and τ fixed. Two things
# follow that the draft quotes. First, with continuation costs paid in steady
# state, profits are lower, so free entry requires a cheaper entry ticket: the
# endogenous economy has disposable firms, the exogenous one expensive ones.
# Second, the exit-flow panel decomposes into a rate and a stock. The exogenous
# specification has a fixed rate, so its flow just tracks the shrinking stock.
# The endogenous one shows a stock falling further yet a flow falling less,
# because the exit rate is climbing as marginal firms cross the threshold.

μ = B.endog.cal.ε / (B.endog.cal.ε - 1)
record!("compB", "BfeEndog", B.endog.ss.ν_f * μ,             "sunk entry cost f_e = ν_f·μ/ρ, endogenous")
record!("compB", "BfeExog",  B.exog.ss.ν_f  * μ,             "same, exogenous")
record!("compB", "BfeRatio", B.exog.ss.ν_f / B.endog.ss.ν_f, "how many times more expensive")
record!("compB", "BdeltaeEndog", 100 * B.endog.ss.δ_e,       "δ_e endogenous, %/month")
record!("compB", "BdeltaeExog",  100 * B.exog.ss.δ_e,        "δ_e exogenous, %/month")

# Exit flow under the (negative) technology shock: impact isolates the χ^c
# margin, because δ_e cannot move in the exogenous specification and N is
# predetermined in both.
ef_endog_z = series(B.endog, :irf_z, :exit_flow; flip=true)
ef_exog_z  = series(B.exog,  :irf_z, :exit_flow; flip=true)
N_endog_z  = series(B.endog, :irf_z, :N;         flip=true)
record!("compB", "BzExitImpactEndog", ef_endog_z[2],          "z: exit flow at h=1, endogenous, %")
record!("compB", "BzExitImpactExog",  ef_exog_z[2],           "z: exit flow at h=1, exogenous, %")
record!("compB", "BzExitPeakEndog",   peak(ef_endog_z).value, "z: exit flow trough, endogenous, %")
record!("compB", "BzExitPeakExog",    peak(ef_exog_z).value,  "z: exit flow trough, exogenous, %")
record!("compB", "BzNPeakEndog",      peak(N_endog_z).value,  "z: N̂ trough, endogenous, %")
# exit_flow = δ̂_e + N̂ by construction in the runners, so the rate contribution
# is the residual. This is the number that makes the decomposition legible.
record!("compB", "BzRateAtPeak",
        ef_endog_z[argmax(abs.(ef_endog_z))] - N_endog_z[argmax(abs.(ef_endog_z))],
        "z: implied δ̂_e at the exit-flow trough, endogenous, %")

# Destruction shock: δ̄ is equated in LEVELS, so the log panel shows the
# exogenous specification jumping twice as far. Recording both makes clear this
# is a normalization, not a behavioral difference.
record!("compB", "BdExitImpactEndog", B.endog.irf_δ.exit_flow[2], "δ: exit flow at h=1, endogenous, %")
record!("compB", "BdExitImpactExog",  B.exog.irf_δ.exit_flow[2],  "δ: exit flow at h=1, exogenous, %")
record!("compB", "BdLevelImpactEndog", B.endog.ss.δ_e * B.endog.irf_δ.exit_flow[2] / 100,
        "δ: level impulse to δ_e, endogenous (should match exogenous)")
record!("compB", "BdLevelImpactExog",  B.exog.ss.δ_e  * B.exog.irf_δ.exit_flow[2]  / 100,
        "δ: level impulse to δ_e, exogenous")
for (tag, case) in (("Endog", B.endog), ("Exog", B.exog))
    record!("compB", "Bd$(tag)Uhl", decay_half_life(series(case, :irf_δ, :u; pp=true)),
            "δ: months for u to halve from peak")
end

# =============================================================================
# 4. Comparison C: variety effects
# =============================================================================
# C sets ρ ≡ 1. The identification is clean because N is predetermined: on
# impact the two specifications must give an identical w^int response, since ρ
# has not yet moved. Everything separating the two paths thereafter is the
# variety externality, with no contamination from the direct effect of z.

w_base_z  = series(C.base,  :irf_z, :w_int; flip=true)
w_novar_z = series(C.novar, :irf_z, :w_int; flip=true)
record!("compC", "CzWintImpactBase",  w_base_z[2],                 "z: w^int at h=1, baseline, %")
record!("compC", "CzWintImpactNovar", w_novar_z[2],                "z: w^int at h=1, no-variety, %")
record!("compC", "CzWintImpactGap",   w_base_z[2] - w_novar_z[2],  "z: the gap at impact (should be 0)")
record!("compC", "CzWintGapYear",
        (w_base_z[11] - w_novar_z[11]) / abs(w_base_z[11]),
        "z: gap as a share of the baseline response at h=10")
record!("compC", "CdWintNovar", peak(C.novar.irf_δ.w_int).value,
        "δ: w^int peak, no-variety (zero by construction, ρ≡1 and z unshocked)")
# ε-driven elasticities in the χ^c decomposition: the market-share effect
# dominates the continuation-value effect whenever ε > 2.
ε = C.base.cal.ε
record!("compC", "ElastProfit", (2 - ε) / (ε - 1), "per-firm profit elasticity w.r.t. N")
record!("compC", "ElastValue",  1 / (ε - 1),       "continuation-value elasticity w.r.t. N")

# =============================================================================
# 5. Comparison D: entry elasticity
# =============================================================================
# ξ plays opposite roles across the two shocks, and what flips the sign is the
# direction in which the shock moves the vacancy value. A technology decline
# lowers Q, so entry falls with the shock and amplifies it. A destruction shock
# shortens a product line's expected life, front-loading the vacancy's value and
# raising Q, so entry leans against the shock and cushions it.

record!("compD", "DxiInvBase", D.base.cal.ξ_inv,        "baseline ξ_inv")
record!("compD", "DxiInvLow",  D.low_ξ.cal.ξ_inv,       "near-free-entry ξ_inv")
record!("compD", "DxiLow",     1 / D.low_ξ.cal.ξ_inv,   "implied ξ = 1/ξ_inv")

Qz_base = peak(series(D.base,  :irf_z, :Q; flip=true)).value
Qz_low  = peak(series(D.low_ξ, :irf_z, :Q; flip=true)).value
ez_base = peak(series(D.base,  :irf_z, :e; flip=true)).value
ez_low  = peak(series(D.low_ξ, :irf_z, :e; flip=true)).value
record!("compD", "DzQBase", Qz_base, "z: Q̂ peak, baseline, %")
record!("compD", "DzQLow",  Qz_low,  "z: Q̂ peak, near-free-entry, %")
record!("compD", "DzERatio", ez_low / ez_base,
        "z: realized amplification of entry (below ξ because Q̂ itself is damped)")

# Destruction shock: the vacancy path is the paper's key robustness result.
# Even at near-free entry, vacancies fall, so the empirically calibrated δ̄_e
# rules out the CK Beveridge-curve shift without needing convex entry costs.
v_low_d  = series(D.low_ξ, :irf_δ, :v; pp=true)
v_base_d = series(D.base,  :irf_δ, :v; pp=true)
record!("compD", "DdVtroughBase", peak(v_base_d).value, "δ: v trough, baseline, pp")
record!("compD", "DdVtroughLow",  peak(v_low_d).value,  "δ: v trough, near-free-entry, pp")
record!("compD", "DdVposMonths",  count(>(0), v_low_d[2:end]) == 0 ? 0 :
        findfirst(<(0), v_low_d[2:end]) - 1,
        "δ: months v stays positive before turning negative, near-free-entry")
record!("compD", "DdUBase", peak(series(D.base,  :irf_δ, :u; pp=true)).value, "δ: u peak, baseline, pp")
record!("compD", "DdULow",  peak(series(D.low_ξ, :irf_δ, :u; pp=true)).value, "δ: u peak, low-ξ, pp")
record!("compD", "DdeltaeTau", D.base.ss.δ_e / TAU_M, "δ_e/τ at the baseline calibration")

# =============================================================================
# Write output
# =============================================================================

open("mechanism_stats.txt", "w") do io
    println(io, "Draft-bound statistics for Section 5.3 and the Comparison D appendix.")
    println(io, "Generated by mechanism_stats.jl from the serialized mechanism IRFs.")
    println(io, "Do not hand-edit; regenerate after re-running the mechanism runners.\n")
    for sec in unique(first.(ROWS))
        println(io, "="^78); println(io, sec); println(io, "="^78)
        for (s, k, v, note) in ROWS
            s == sec && @printf(io, "  %-22s %14.5f   %s\n", k, v, note)
        end
        println(io)
    end
end

# LaTeX macros, so the draft can \input this file and stop carrying hard-coded
# numbers. Macro names are letters-only as TeX requires. Trailing comments are
# stripped to ASCII: this file gets \input into the manuscript, and stray
# Unicode in a comment is a needless way to break someone else's build.
ascii_note(s) = replace(s,
    "δ" => "delta", "ω" => "omega", "τ" => "tau", "ξ" => "xi", "χ" => "chi",
    "ε" => "eps", "μ" => "mu", "ρ" => "rho", "θ" => "theta", "ψ" => "psi",
    "ν" => "nu", "σ" => "sigma", "̂" => "", "^" => "", "≡" => "=", "→" => "->",
    "·" => "*", "≤" => "<=", "≥" => ">=")

open("mechanism_stats.tex", "w") do io
    println(io, "% Auto-generated by mechanism_stats.jl. Do not edit by hand.")
    println(io, "% \\input{mechanism_stats} in the preamble, then use e.g. \\SSu.")
    for (_, k, v, note) in ROWS
        @printf(io, "\\newcommand{\\%s}{%.4g}%% %s\n", k, v, ascii_note(note))
    end
end

println("Wrote mechanism_stats.txt and mechanism_stats.tex ($(length(ROWS)) quantities)")
