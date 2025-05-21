include("steady_state.jl")

targets = (labor_share=0.66, 
           dest_ann=0.10, 
           r_ann=0.04, 
           f =0.41, 
           η_L=0.6, #elast. of matching function
           q=0.8, 
           sep=0.031, 
           b_ratio=0.71, 
            x_v=0.20, 
            ξ_inv=1, # congestion elasticity
            ε=4.3, # elasticity of sub.
            σ=1.0, # log utility
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
            "Vacancy share",
            "Investment in new product lines",
            "Sunk vacancy cost share",
            "Entrant share",
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

output = IOBuffer()
show(output, MIME("text/latex"),df_shares)
df_shares_table = String(take!(output))
print(df_shares_table)