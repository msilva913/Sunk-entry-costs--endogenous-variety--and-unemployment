# cutoff_profit_share.jl
# =============================================================================
# rem:dilution_variety. The survival cutoff (eq:cutoff_eq) is
#   χ^c = (μ−1)Y^c/(μN) + f_e ρ(N)/μ        (flow profit + value of a line).
# Holding zL^c fixed, flow profit has elasticity 1/(ε−1) − 1 in N (variety and dilution)
# and the value of a line 1/(ε−1) (variety only), so the cutoff has elasticity
#   1/(ε−1) − s_R,   s_R = flow-profit share of the cutoff.
# It is positive (a smaller mass of lines lowers the cutoff and induces exit) for
# ε < 1 + 1/s_R. Reports s_R at the code defaults (D1, PATH B) for several ε.
# Output: cutoff_profit_share.txt
# =============================================================================
include("steady_state.jl")
using Printf

lines = String[]
log!(s) = (push!(lines, s); println(s))
log!("Flow-profit share s_R of the survival cutoff at the code defaults")
for ε in (4.3, 3.0, 2.4, 2.0)
    cal = calibrate_shares(merge(TARGETS, (ε = ε,)))
    st  = steady_state(cal)
    μ   = ε / (ε - 1)
    R   = (μ - 1) * st.Y_c / (μ * st.N)
    V   = cal.f_e * st.ρ / μ
    @assert abs(R + V - st.x_c) < 1e-8 * st.x_c      # eq:cutoff_eq holds in steady state
    sR  = R / (R + V)
    log!(@sprintf("  ε = %.1f: flow profit %.4f, value of a line %.3f, cutoff %.3f; s_R = %.4f; elasticity 1/(ε−1) − s_R = %.4f; positive for ε < 1 + 1/s_R = %.1f",
                  ε, R, V, st.x_c, sR, 1 / (ε - 1) - sR, 1 + 1 / sR))
end
write("cutoff_profit_share.txt", join(lines, "\n") * "\n")
