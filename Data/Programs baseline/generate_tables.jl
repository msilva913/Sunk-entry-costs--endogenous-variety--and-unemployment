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
    @unpack f_e, δ, s, z, b, ϕ, σ, ε, A, η_L, κ, ξ_inv, f_m, ψ, p_0 = cal
    @unpack X_Y, Xc_Y, dest_ann, dest_end_frac, f, q, sep, x_v, r_ann, N, w = targets

    r = (1 + r_ann)^(1/12) - 1
    β = 1 / (1 + r)
    δ_e = 1 - (1 - dest_ann)^(1/12)

    Parameter = [
        L"b", L"\kappa", L"\xi^{-1}", L"\sigma", L"\delta_{\text{end,share}}", L"p_0",  # Estimated
        L"r", L"\eta_L", L"\varepsilon", L"\delta_e",                        # Directly set
        L"z", L"f_e",                                                        # Dependent: normalizations
        L"s", L"\phi", L"A", L"f_m", L"\psi"                      # Dependent: long-run targets
    ]

    Targets = [
        "Estimated", "Estimated", "Estimated", "Estimated", "Estimated", "Estimated",
        "Real interest rate (annual)", "Elasticity of matching function", "Elasticity of substitution", "Establishment exit rate (annual)",
        "Steady-state wage", "Steady-state mass of firms",
        "Aggregate separation rate",
        "Recruiting cost share: " * L"$X/Y$",
        "Job finding rate",
        L"Endogenous destruction share: " * L"$(1-F(x^c))/\delta_e$",
        "Fixed cost share: " * L"$X^c/Y$"
    ]

    Value = [
        "-", "-", "-", "-", "-", "-",
        r_ann, η_L, ε, dest_ann,
        w, N,
        sep, X_Y, f, dest_end_frac, Xc_Y
    ]

    Calibration = [
        b, κ, ξ_inv, σ, dest_end_frac, p_0,
        r, η_L, ε, δ_e,
        z, f_e,
        s, ϕ, A, f_m, ψ
    ]

    df = DataFrame(
        Parameter = Parameter,
        Targets = Targets,
        Value = map(x -> isa(x, Number) ? round(x, sigdigits = 2) : x, Value),
        Calibration = round.(Calibration, sigdigits = 3)
    )
    return df
end


targets = (
    X_Y=0.015,        # recruiting cost share of output
    Xc_Y=0.10,         # fixed cost share of output (Abraham, Bormans, Konings, Roeger)
    dest_ann=0.0754,   # annual product destruction rate
    dest_end_frac=0.5, # endogenous share of destruction rate (Estimated)
    p_0=0.5,           # probability of drawing from
    #dest_el=1.0,      # Destruction elasticity wrt x_c
    f=0.41,            # job-finding rate 
    η_L=0.6,           # elasticity of matching fun wrt unemployment
    q=0.8,             # vacancy filling rate,
    sep=0.031,         # aggregate separation rate
    b_ratio=0.71,      # ratio of unemployment benefits to wage (Estimated)
    x_v=1.0,           # (Estimated)
    ξ_inv=1,           # entry elasticity inverse (Estimated)
    r_ann=0.04,        # annual discount rate
    ε=4.3,             # elasticity of substitution (BGM, Compustat)
    σ=1.0,             # inverse IES (Estimated)
    N=1.0,             # SS mass of firms (normalization)
    w=1.0              # SS wage (normalization)
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

# Text that can be exported to WhatsApp
show(output, MIME"text/plain"(), df; allrows=true, allcols=true)
text = "" * String(take!(output)) * "" 
clipboard(text) # copy to system clipboard
println("Copied ", length(text), " chars.")

function shares_table(ss::NamedTuple)
    # Define additional variables 
    @unpack w_int, w, L, X, Y, N, d_f, δ_e, dest_end_frac, entrant_vac_share = ss
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
         ss.X_c/ss.Y, ss.dest_end_frac,  ss.entrant_vac_share, ss.labor_share, recruiter_profit_share, retailer_profit_share,
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