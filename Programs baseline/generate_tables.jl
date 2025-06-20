"""
Generate tables for 
    1) calibration 
    2) steady-state shares
"""

include("steady_state.jl")
cd("C:/Users/msilva913/Documents/GitHub/Sunk_entry_costs_endogenous_variety_unemployment")
# Load posterior mode 
posterior_mode = matopen("posterior_mode.mat")
posterior_mode = read(posterior_mode, "posterior_mode")

targets = (labor_share=0.66, 
           dest_ann=0.10, 
           r_ann=0.04, 
           f =0.41, 
           η_L=0.6, #elast. of matching function
           q=0.8, 
           sep=0.031, 
           b_ratio=posterior_mode["b_ratio"], 
            x_v=posterior_mode["x_v"], 
            ξ_inv=posterior_mode["xi_inv"], # congestion elasticity
            ε=posterior_mode["epsi"], # elasticity of sub.
            σ=posterior_mode["sigma"], # log utility
            N=1, w=1.0)

cal = calibrate_labor_share(targets)
ss = steady_state(cal)
df = calibration_table(cal, targets)

# Generate latex output
output = IOBuffer()
show(output, MIME("text/latex"),df)
df_table = String(take!(output))
print(df_table)

function shares_table(ss::NamedTuple)
    df = DataFrame(
        Share = [
            "Consumption share",
            "Recruiting cost share",
            "Investment in new product lines",
            "Sunk vacancy cost share",
            "New vacancy share",
            "Search wedge",
            "Market power wedge",
            "Stock market cap to GDP"
        ],
        Symbol = [L"C/Y", L"X/Y", L"\nu N_e/Y", L"X_v/Y", L"e/v", L"w/w^R", L"w^RL/Y", L"M/(12*Y)"],
        Value = round.([ss.cons_share, ss.vacancy_share, ss.inv_new_firm_share, ss.sunk_vac_cost_share, ss.entrant_share, ss.search_wedge, ss.recruiter_share, ss.M/(12*ss.Y)], sigdigits=2),

    )
    return df
end

df_shares = shares_table(ss)

"""
output = IOBuffer()
show(output, MIME("text/latex"),df_shares)
df_shares_table = String(take!(output))
print(df_shares_table)
"""
function dataframe_to_latex(df::DataFrame)
    header = [
        "\\begin{table}[ht]",
        "\\centering",
        "\\begin{tabular}{l|cc}",
        "\\toprule",
        "Share & Symbol & Value \\\\ \\midrule"
    ]

    rows = [
        string(
            row.Share, " & ", row.Symbol, " & ", round(row.Value, digits=3), " \\\\"
        ) for row in eachrow(df)
    ]

    footer = [
        "\\bottomrule",
        "\\end{tabular}",
        "\\caption{}",
        "\\end{table}"
    ]

    return join([header; rows; footer], "\n")
end

latex_table = dataframe_to_latex(df_shares)
println(latex_table)