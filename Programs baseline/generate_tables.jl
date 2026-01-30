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
    @unpack f_e, δ, s, z, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, f_m, ψ, s = cal
    @unpack X_Y, Xc_Y, dest_ann, dest_end_frac, f, η_L, q, sep, b_ratio, x_v, ξ_inv, ε, r_ann, σ, N, w = targets

    r = (1+r_ann)^(1/12)-1
    β = 1/(1+r)
    δ_e = 1-(1-dest_ann)^(1/12)
    μ = ε/(ε-1)
    τ = (1-δ_e)*(1-s)
    # Creating a DataFrame for the table
    df = DataFrame(
        Parameter = [L"r", L"\eta_L", L"b", L"\delta", L"\xi^{-1}", L"\varepsilon", L"\sigma", L"\kappa", L"z", L"f_e", L"s",
         L"\phi", L"A", L"f_m", L"\psi"],
        Targets = [
            "Real interest rate",
            "Elasticity of matching function",
          "Estimated",   #  "Replacement ratio b/w",
            "Job destruction from obsolescence",
           "Estimated", # "Elasticity of vacancy value",
            "Markup", # Elasticity of substitution
          # "Investment share",
           "Estimated",  # "Risk aversion",
           "Estimated", # "Share of sunk vacancy costs to overall hiring costs",
            "Steady-state wage",
            "Steady-state mass of firms",
            "Match separation rate",
            "Recruiting cost share: X/Y",
            "Job finding rate",
            L"Endogenous destruction share: "*L"$(1-F(x_c))/\delta_e$",
            "Fixed cost share: "*L"$X^c/Y$",
        ],
        Value = map(x -> isa(x, Number) ? round(x, sigdigits=2) : x, [r_ann, η_L, "-", δ, "-", ε, "-", "-", w, N, τ, X_Y, 0.41, 0.5, 0.20]),
        Calibration = round.([r, η_L, b, δ, ξ_inv, ε, σ, κ, z, f_e, s, ϕ, A, f_m, ψ], sigdigits=3)
    )

    # Save the DataFrame as a PDF table
    return df
end

targets = (
    X_Y=0.015,        # recruiting cost share of output
    Xc_Y=0.20,         # fixed cost share of output (Abraham, Bormans, Konings, Roeger)
    dest_ann=0.0754,   # annual product destruction rate
    dest_end_frac=0.5, # endogenous share of destruction rate
    f=0.41,            # job-finding rate, 
    η_L=0.6,           # elasticity of matching fun wrt unemployment
    q=0.8,             # vacancy filling rate,
    sep=0.031,         # aggregate separation rate , 
    b_ratio=0.71,      # ratio of unemployment benefits to wage,
    x_v=1.0, 
    ξ_inv=1, 
    r_ann=0.04,        # annual discount rate
    ε=4.3,             # Elasticity of substitution (BGM, Compustat)
    σ=1.0,             # Inverse IES
    N=1.0,             # SS mass of forms (normalization)
    w=1.0,             # SS wage (normalization)
    #ψ=1.5)
)

"""
targets = (X_Y = 0.015, # vacancy share target,  influences ϕ
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
"""

cal = calibrate_shares(targets)
ss = steady_state(cal)

df = calibration_table(cal, targets)

# Generate latex output
output = IOBuffer()
show(output, MIME("text/latex"),df)
df_table = String(take!(output))
print(df_table)

function shares_table(ss::NamedTuple)
    # Define additional variables 
    @unpack w_int, w, L, X, Y, N, d_f, δ_e, end_dest_share, entrant_vac_share = ss
    recruiter_profit_share = ((w_int-w)*L-X)/Y
    retailer_profit_share =  (N*d_f)/Y

    df = DataFrame(
        Share = [
            "Annual interest rate",
            "Gross markup",
            "Consumption share",
            "Vacancy rate",
            "Unemployment rate",
            "Market tightness",
             "Business formation share",
            "Recruiting cost share",
            "Fixed cost share",
            "Endogenous share of destruction shock",
            "New vacancy share",
            #"Search wedge",
            #"Market power wedge",
            "Labor share",
            "Recruiter profit share",
            "Retailer profit share",
            "Value of a vacancy",
            "Value of a filled job",
            "Stock market cap to GDP"
        ],
        Symbol = [L"(1+r)^12-1", L"\mu", L"C/Y", L"v", L"u", L"\theta",  L"\nu N_e/Y", L"X/Y", L"Xc/Y", L"(1-F(x^c))/\delta_e",
       L"e/v", L"wL/Y", L"(w^{int}-w)*L-X)/Y", L"N*d_f/Y",  L"Q", L"J", L"M/(12*Y)"],
        Value = round.([ss.ann_int_rate, ss.μ, ss.cons_share, ss.v, ss.u, ss.θ, ss.inv_new_firm_share, ss.vacancy_share,
         ss.X_c/ss.Y, ss.end_dest_share,  ss.entrant_vac_share, ss.labor_share, recruiter_profit_share, retailer_profit_share,
        ss.Q, ss.J, ss.M/(12*ss.Y)], sigdigits=3),

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