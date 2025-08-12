"""
Generate tables for 
    1) calibration 
    2) steady-state shares
"""

include("steady_state.jl")
using MAT
#cd("C:/Users/msilva913/Documents/GitHub/Sunk_entry_costs_endogenous_variety_unemployment")
# Load posterior mode 
#posterior_mode = matopen("posterior_mode.mat")
#posterior_mode = read(posterior_mode, "posterior_mode")

function calibration_table(cal, targets)
    @unpack f_e, τ, δ, z, b, ϕ, ρ, σ, ε, A, η_L, ξ_inv, A, F, x_m, κ, s = cal
    @unpack labor_share, dest_ann, f, η_L, q, sep, b_ratio, x_v, ξ_inv, ε, r_ann, σ, N, w = targets

    #ρ = (1+r_ann)^(1/12)-1
    β = 1/(1+ρ)
    δ = 1-(1-dest_ann)^(1/12)
    μ = ε/(ε-1)

    # Creating a DataFrame for the table
    df = DataFrame(
        Parameter = [L"\rho", L"\eta_L", L"b", L"\delta", L"\xi^{-1}", L"\varepsilon", L"\sigma", L"\kappa", L"z", L"f_e", L"s",
         L"\phi", L"A", L"F"],
        Targets = [
            "Real interest rate",
            #"Recruiting cost share",
            "Elasticity of matching function",
          #  "Replacement ratio b/w",
          "Estimated",
            "Job destruction from obsolescence",
           # "Elasticity of vacancy value",
           "Estimated",
           # "Markup",
           "Investment share",
           # "Risk aversion",
           "Estimated",
           # "Share of sunk vacancy costs to overall hiring costs",
           "Estimated",
            "Steady-state wage",
            "Steady-state mass of firms",
            "Aggregate separation rate",
            "Labor share",
            "Job finding rate",
            "Vacancy filling rate"
        ],
        Value = map(x -> isa(x, Number) ? round(x, sigdigits=2) : x, [r_ann, η_L, "-", δ, "-", ε, "-", "-", w, N, τ, labor_share, 0.41, 0.80]),
        Calibration = round.([ρ, η_L, b, δ, ξ_inv, ε, σ, κ, z, f_e, s, ϕ, A, F], sigdigits=3)
    )

    # Save the DataFrame as a PDF table
    return df
end

"""
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
"""
targets = (labor_share=0.66, # influences ϕ
           dest_ann=0.0754, #21% of job destruction from obsolescence 
           f =0.41, # fixed, turnover means
           η_L=0.6, # based on time-series regressions
           q=0.8, # fixed, turnover means
           sep=0.031, # imputed from unemployment flows according to Shimer (2005)
           b_ratio=0.71, #estimated
           x_v=1.0, # estimated, affects X/Y, 
           ξ_inv=1.0, # estimated, affects X/Y, 
           #ξ_inv = 1/0.265,
           #X_Y=0.015, #vacancy share target, as Shao and Silos
           r_ann=0.04, #4% annual interest rate
           #C_Y=0.80, # consumption share, influences value of ε
           ε = 4.3,
           σ=1.0, # benchmark corresponding to log preferences, estimated
           N=1.0, # normalization: pins down f_e
           w=1.0, # normalization: we express values relative to wage
)

cal = calibrate_shares(targets)
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
            "Annual interest rate",
            "Gross markup",
            "Consumption share",
            "Vacancy rate",
            "Unemployment rate",
            "Market tightness",
            "Recruiting cost share",
            "Business formation share",
            "Sunk vacancy cost share",
            "New vacancy share",
            "Search wedge",
            "Market power wedge",
            "Value of a vacancy",
            "Value of a filled job",
            "Stock market cap to GDP"
        ],
        Symbol = [L"(1+ρ)^12-1", L"\mu", L"C/Y", L"v", L"u", L"\theta", L"X/Y", L"\nu N_e/Y", L"X_v/Y", L"e/v", L"w/w^R", L"w^RL/Y", 
        L"Q", L"J", L"M/(12*Y)"],
        Value = round.([ss.ann_int_rate, ss.μ, ss.cons_share, ss.v, ss.u, ss.θ, ss.vacancy_share, ss.inv_new_firm_share, 
        ss.sunk_vac_cost_share, ss.entrant_share, ss.search_wedge, ss.recruiter_share, ss.Q, ss.J, ss.M/(12*ss.Y)], sigdigits=3),

    )
    return df
end

df_shares = shares_table(ss)
pprint(df_shares)
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