# ags_xi_sweep.jl
# =============================================================================
# Exploratory, not draft-bound. The trade-off in the elasticity of vacancy creation. Elastic
# creation (small ξ_inv) amplifies z (ags_kappa_sweep.jl), but some inelasticity is needed for
#   (1) a mild comovement of unemployment with labor productivity under z shocks (the point in
#       Gabrovski and Silva): with inelastic creation, u lags productivity, which weakens the
#       contemporaneous correlation;
#   (2) Beveridge-curve dynamics after δ shocks: destroyed vacancies are replaced slowly, so v
#       falls while u rises.
# Grid: ξ_inv ∈ XI_GRID at x_v = 0.5 (the default) and 0.25 (κ's share of a hire 75 %). Everything
# else at AGS_TARGETS (Λ̄_r = 0.25, p_r = 0.8, ψ_r = 1, X/Y 1.5 %).
#
# Moments come from time_series_fun.jl moment_table (the form of the draft's Table 3), on
# quarterly averages of the simulated monthly series, HP-1600, with all shocks and with each shock
# alone. Model counterparts of the data series: ags_model.jl model_quarterly.
# Data: data/results/q3_data_svar_cycles_hp.csv from q3_data_svar.py (same HP cycles as the SVAR),
# 1951Q1–2019Q4 and SVAR A's window 1992Q3–2019Q4. δ starts in 1992Q3 (pairwise complete).
#
# Output: ags_xi_sweep.txt.
# =============================================================================

include("ags_model.jl")
using CSV

const XI_GRID  = (4.0, 2.0, 1.0, 0.5, 0.25, 0.1, 0.03, 0.01)
const XV_CASES = (0.5, 0.25)
const FULL_AT  = ((1.0, 0.5), (0.01, 0.5))        # (ξ_inv, x_v) cases with full tables
const TABLE_VARS = ["u", "v", "s", "f", "δ", "z"]  # Table 3 order (τ there is s here)
const LAGS = -8:8                                  # quarters, for cor(u_t, z_{t−k})
const DATA_CSV = joinpath(@__DIR__, "..", "Bartek analysis", "data", "results",
                          "q3_data_svar_cycles_hp.csv")

out = String[]
note(str) = (push!(out, str); println(str))
"""Simulated quarterly cycles with the shocks in `cols` (1 z, 2 δ, 3 s)."""
function simulate_cycles(r, cols)
    eta_m = hcat([shock_col(i, SDS[i]) for i in cols]...)
    sim = DataFrame(simulate_model((; model..., ne = length(cols)), r.sol_mat, T_SM, eta_m,
                                   r.SS, false, true), varnames)
    return model_cycles(sim, r.cal)
end

"""Lag k (quarters) maximizing |cor(u_t, z_{t−k})|, and that correlation; k > 0: u lags z."""
function peak_lag(u, z)
    best_k, best_c = 0, 0.0
    for k in LAGS
        if k >= 0
            c = cor(u[1+k:end], z[1:end-k])
        else
            c = cor(u[1:end+k], z[1-k:end])
        end
        if abs(c) > abs(best_c)
            best_k, best_c = k, c
        end
    end
    return best_k, best_c
end

"""
Series a one-shock simulation actually moves. A shock that does not reach a series leaves it at
rounding noise (σ ≈ 1e-17), whose correlations are meaningless, so it is dropped from the table.
"""
active(c) = [v for v in TABLE_VARS if std(c[!, v]) > 1e-8]

"""cor(a, b) from a moment_table, whichever triangle holds it (the table fills only one)."""
function table_cell(tab, a, b)
    x = tab[findfirst(==(a), tab.variable), b]
    if ismissing(x)
        x = tab[findfirst(==(b), tab.variable), a]
    end
    return x
end

# ── data ────────────────────────────────────────────────────────────────────
raw = CSV.read(DATA_CSV, DataFrame; missingstring = "")
rename!(raw, :jf => :f, :delta => :δ, :theta => :θ)
note("AGS ξ_inv sweep (ags_xi_sweep.jl). Quarterly HP-1600 cycles; moments in the form of Table 3.")
data_rows = Dict{String,Any}()
for (lab, sel) in (("1951Q1–2019Q4", trues(nrow(raw))), ("1992Q3–2019Q4 (SVAR A)", raw.quarter .>= "1992Q3"))
    d = raw[sel, :]
    tab = moment_table(d[!, TABLE_VARS]; rel_to = :z)
    full = .!ismissing.(d.u) .& .!ismissing.(d.z)
    k, c = peak_lag(Float64.(d.u[full]), Float64.(d.z[full]))
    θz = std(skipmissing(d.θ)) / std(skipmissing(d.z))
    note("\nData, $lab")
    foreach(note, moment_table_text(tab))
    note(@sprintf("  σ(θ)/σ(z) %.2f; peak |cor(u_t, z_{t−k})| at k = %d: %+.3f", θz, k, c))
    data_rows[lab] = (tab = tab, k = k, c = c)
end

# ── sweep ───────────────────────────────────────────────────────────────────
note("\nSweep. 'all' = all shocks; '|z', '|δ', '|s' = that shock alone. lag k: quarters by which u")
note("lags z at the peak |cor(u_t, z_{t−k})| under z alone. δ IRF: u peak, v trough (%, month).")
full_tables = String[]
for xv in XV_CASES
    note(@sprintf("\nx_v = %.2f (κ's share of the cost of a hire %d %%)", xv, round(Int, 100 * (1 - xv))))
    note(@sprintf("  %6s %6s | %5s %5s %5s %5s | %6s %6s %3s %6s | %6s %6s %6s %6s | %6s | %13s %13s",
                  "ξ_inv", "ϕ", "σu z", "σu δ", "σu s", "all", "cor uz", "|z", "k", "at k",
                  "cor uv", "|z", "|δ", "|s", "σθ/σz", "δ: u peak", "δ: v trough"))
    for ξ in XI_GRID
        local r = nothing
        try
            r = solve_ags(merge(AGS_TARGETS, (x_v = xv, ξ_inv = ξ)))
        catch err
            note(@sprintf("  %6.2f failed: %s", ξ, first(sprint(showerror, err), 90)))
            continue
        end
        c_all = simulate_cycles(r, 1:3)
        c_z, c_d, c_s = simulate_cycles(r, [1]), simulate_cycles(r, [2]), simulate_cycles(r, [3])
        t_all = moment_table(c_all[!, TABLE_VARS]; rel_to = :z)
        t_z   = moment_table(c_z[!, active(c_z)]; rel_to = :z)
        t_d   = moment_table(c_d[!, active(c_d)])
        t_s   = moment_table(c_s[!, active(c_s)])
        k, ck = peak_lag(c_z.u, c_z.z)
        ird   = irf(r, 2)
        pu, pv = ird[2:end, :u], ird[2:end, :v]
        iu, iv = argmax(pu), argmin(pv)
        note(@sprintf("  %6.2f %6.3f | %5.2f %5.2f %5.2f %5.2f | %+6.3f %+6.3f %3d %+6.3f | %+6.3f %+6.3f %+6.3f %+6.3f | %6.2f | %6.3f (%3d) %6.3f (%3d)",
                      ξ, r.cal.ϕ, 100std(c_z.u), 100std(c_d.u), 100std(c_s.u), 100std(c_all.u),
                      table_cell(t_all, "u", "z"), table_cell(t_z, "u", "z"), k, ck,
                      table_cell(t_all, "v", "u"), table_cell(t_z, "v", "u"),
                      table_cell(t_d, "v", "u"), table_cell(t_s, "v", "u"),
                      std(c_all.θ) / std(c_all.z), pu[iu], iu - 1, pv[iv], iv - 1))
        if (ξ, xv) in FULL_AT
            push!(full_tables, @sprintf("\nModel, ξ_inv = %.2f, x_v = %.2f", ξ, xv))
            for (lab, t) in (("all shocks", t_all), ("z alone", t_z), ("δ alone", t_d), ("s alone", t_s))
                push!(full_tables, "  $lab")
                append!(full_tables, moment_table_text(t))
            end
        end
    end
end
note("\nσu: σ(u) in %. Data (SVAR A window, ad hoc from shares): σu z 6.32, δ 3.64, s 3.82; three 8.24.")
note("\n" * "-"^78 * "\nFull moment tables (Table 3 form)")
foreach(note, full_tables)

open("ags_xi_sweep.txt", "w") do io
    foreach(l -> println(io, l), out)
end
