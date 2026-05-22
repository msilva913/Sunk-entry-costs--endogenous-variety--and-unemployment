# export_endog_exit_irfs.jl
# =============================================================================
# Load irf_endog_exit.jls and write all IRFs + SS values to CSV.
# Run from the Programs baseline/ directory:
#   julia export_endog_exit_irfs.jl
# Produces four files:
#   irf_endog_z.csv, irf_endog_delta.csv
#   irf_exog_z.csv,  irf_exog_delta.csv
# and a summary table:   ss_comparison.csv
# =============================================================================

using Serialization, DataFrames, CSV

output_B = deserialize("irf_endog_exit.jls")

# ── IRF DataFrames ────────────────────────────────────────────────────────────
CSV.write("irf_endog_z.csv",     output_B.endog.irf_z)
CSV.write("irf_endog_delta.csv", output_B.endog.irf_δ)
CSV.write("irf_exog_z.csv",      output_B.exog.irf_z)
CSV.write("irf_exog_delta.csv",  output_B.exog.irf_δ)
println("Wrote: irf_endog_z.csv, irf_endog_delta.csv, irf_exog_z.csv, irf_exog_delta.csv")

# ── Steady-state comparison table ─────────────────────────────────────────────
endog_ss = output_B.endog.ss
exog_ss  = output_B.exog.ss
endog_cal = output_B.endog.cal
exog_cal  = output_B.exog.cal

ss_fields = [:u, :v, :θ, :δ_e, :N, :N_e, :x_c, :K, :Q, :e, :w_int,
             :labor_prod, :π_s, :dest_el, :dest_end_frac]

rows = []
for f in ss_fields
    e_val = hasproperty(endog_ss, f) ? getfield(endog_ss, f) : missing
    x_val = hasproperty(exog_ss,  f) ? getfield(exog_ss,  f) : missing
    push!(rows, (variable=string(f), endog=e_val, exog=x_val))
end
# Add calibrated parameters not in ss
for (f, getter) in [(:s, c->c.s), (:κ, c->c.κ), (:ϕ, c->c.ϕ),
                    (:ψ, c->c.ψ), (:f_e, c->c.f_e), (:z, c->c.z)]
    push!(rows, (variable=string(f),
                 endog=getter(endog_cal),
                 exog=getter(exog_cal)))
end

ss_df = DataFrame(rows)
CSV.write("ss_comparison.csv", ss_df)
println("Wrote: ss_comparison.csv")

# ── Quick peak-response table ──────────────────────────────────────────────────
println("\n── Peak responses (100×log dev unless noted) ───────────────────")
for (shock, irf_e, irf_x) in [("z",  output_B.endog.irf_z, output_B.exog.irf_z),
                               ("δ",  output_B.endog.irf_δ, output_B.exog.irf_δ)]
    println("  $shock shock:")
    for col in [:u, :v, :θ, :N, :N_e, :K, :e, :w_int, :exit_flow]
        if hasproperty(irf_e, col) && hasproperty(irf_x, col)
            pe = round(maximum(abs.(irf_e[!, col])), digits=4)
            px = round(maximum(abs.(irf_x[!, col])), digits=4)
            println("    $(rpad(string(col),12)) endog=$(pe)   exog=$(px)")
        end
    end
end
