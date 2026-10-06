# predetermined_vacancy_share.jl
# =============================================================================
# Computes the share of the steady-state vacancy stock that is PREDETERMINED
# (inherited from t-1 and therefore exposed to a destruction shock), and the
# entry elasticity a destruction shock must induce for total vacancies to RISE.
#
# WHY
# ---
# The log-linear vacancy law of motion (draft eq:v_ll) gives, on impact of a
# pure delta shock with v_{t-1} and s at steady state,
#
#     dv_t = de_t - v_pre * ddelta_e,t ,      v_pre = v - e
#
# so vacancies rise iff new postings exceed the destroyed predetermined stock.
# ⚠ Oct 6, 2026: the next display is the bar per UNIT LEVEL of delta_e, not per 1%.
# Per 1% (both log deviations) it is e-hat/delta-hat > delta_e*v_pre/e = v_pre/(v+L).
# See beveridge_free_entry.jl and the rewritten note. Original text follows.
# In log deviations (e-hat = de/e, delta-hat = ddelta_e/delta_e) that condition is
#
#     e-hat / delta-hat > v_pre / e  ==  v/e - 1 .
#
# Using the steady-state entry relation e = delta_e*(v + L) (eq:e_ss,
# steady_state.jl:217) this threshold is
#
#     T(delta_e) = [ v/(v+L) ] / delta_e - 1 ,
#
# i.e. inversely proportional to the destruction rate. Equivalently, since the
# f, q and sep targets pin v/(v+L) and tau,
#
#     T = kappa_v / (delta_e/tau) - 1 ,   kappa_v = [v/(v+L)] / tau ,
#
# so the threshold is governed entirely by delta_e/tau. Coles-Kelishomi set
# delta_e = tau, the MAXIMUM admissible value (s >= 0 requires delta_e <= tau),
# which minimizes the threshold over the whole admissible range.
#
# SELF-CONTAINED. Replicates the relevant steady-state algebra of
# steady_state.jl without the nonlinear solve, since every quantity here is
# closed-form in the targets. Target values and the two helper identities are
# mirrored from steady_state.jl and asserted against it where possible:
#   compute_separation_rate(d,s) = 1 - (1-d)*(1-s)      [steady_state.jl:41]
#   compute_unemployment(t,f,d)  = t / (t + (1-d)*f)    [steady_state.jl:65]
#   theta = f_corr/q_corr = f/q                        [Stage 1]
#   e = delta_e*(v + L)                                [steady_state.jl:217]
#
# OUTPUT
#   predetermined_vacancy_share.txt
# =============================================================================

using Printf

# --- Targets, mirrored from steady_state.jl TARGETS --------------------------
const F_GROSS = 0.41    # gross job-finding rate [JOLTS]
const Q_GROSS = 0.80    # gross vacancy-filling rate [JOLTS]
const SEP     = 0.031   # aggregate separation rate [Shimer 2005 / JOLTS]

# Destruction calibrations to compare (monthly delta_e implied by annual rate)
delta_e_from_ann(dest_ann) = 1 - (1 - dest_ann)^(1 / 12)

struct Calib
    name::String
    delta_e::Float64
    note::String
end

const CALIBS = [
    Calib("D1 settled (BED Deaths, 3.20%/yr)", delta_e_from_ann(0.0320),
          "decisions.md D1; delta_e/tau = 0.087"),
    Calib("Code baseline (superseded, 7.54%/yr)", delta_e_from_ann(0.0754),
          "steady_state.jl TARGETS as of Oct 2026; delta_e/tau = 0.21"),
    Calib("BGM / GS / Shao-Silos (10%/yr)", delta_e_from_ann(0.10),
          "literature comparison, see Notes/delta_calibration..."),
    Calib("Coles-Kelishomi (delta_e = tau)", SEP,
          "all separation is destruction; s = 0; MAXIMUM admissible delta_e"),
]

# --- Steady state, closed form ----------------------------------------------
"""
    ss_block(delta_e)

Returns the steady-state objects that the predetermined-share calculation
needs. theta, u, v and L are invariant to delta_e, because the f, q and sep
targets are held fixed and s adjusts to absorb the change (this is the design
of Comparison A: "all other targets unchanged -> same SS u, v, theta").
"""
function ss_block(delta_e)
    tau = SEP
    s   = (tau - delta_e) / (1 - delta_e)          # inverts eq. at :41
    f_c = F_GROSS / (1 - delta_e)                  # Stage 1 gross correction
    q_c = Q_GROSS / (1 - delta_e)
    theta = f_c / q_c                              # = F_GROSS/Q_GROSS
    u   = tau / (tau + (1 - delta_e) * f_c)        # = tau/(tau+F_GROSS)
    L   = 1 - u
    v   = theta * u
    e   = delta_e * (v + L)                        # eq:e_ss
    v_pre = v - e
    return (; tau, s, theta, u, L, v, e, v_pre)
end

lines = String[]
push!(lines, "="^104)
push!(lines, "PREDETERMINED VACANCY SHARE AND THE ENTRY THRESHOLD FOR A DESTRUCTION SHOCK")
push!(lines, "Targets held fixed: f = $(F_GROSS), q = $(Q_GROSS), sep = $(SEP) (monthly)")
push!(lines, "="^104)
push!(lines, "")

b0 = ss_block(delta_e_from_ann(0.0320))
kappa_v = (b0.v / (b0.v + b0.L)) / b0.tau
push!(lines, @sprintf("theta = %.6f   u = %.6f   L = %.6f   v = %.6f   (invariant to delta_e)",
                      b0.theta, b0.u, b0.L, b0.v))
push!(lines, @sprintf("vacancy share of all positions  v/(v+L) = %.6f", b0.v / (b0.v + b0.L)))
push!(lines, @sprintf("kappa_v = [v/(v+L)]/tau         = %.6f", kappa_v))
push!(lines, "")
push!(lines, "Threshold entry elasticity for total vacancies to RISE on impact:")
push!(lines, "    T = v_pre/e = v/e - 1 = [v/(v+L)]/delta_e - 1 = kappa_v/(delta_e/tau) - 1")
push!(lines, "")

hdr = @sprintf("%-40s %10s %9s %9s %9s %9s", "Calibration", "delta_e", "d_e/tau",
               "e/v", "v_pre/v", "T")
push!(lines, hdr)
push!(lines, "-"^104)

for c in CALIBS
    b = ss_block(c.delta_e)
    dtau = c.delta_e / b.tau
    ev   = b.e / b.v
    vpv  = b.v_pre / b.v
    T    = b.v_pre / b.e
    # Cross-check the two closed forms agree
    @assert isapprox(T, (b.v / (b.v + b.L)) / c.delta_e - 1; rtol = 1e-10)
    @assert isapprox(T, kappa_v / dtau - 1; rtol = 1e-10)
    @assert isapprox(vpv, 1 - dtau / kappa_v; rtol = 1e-10)
    push!(lines, @sprintf("%-40s %10.6f %9.4f %9.4f %9.4f %9.3f",
                          c.name, c.delta_e, dtau, ev, vpv, T))
end
push!(lines, "-"^104)
push!(lines, "")
push!(lines, "Columns: delta_e monthly destruction rate; d_e/tau its share of total separations;")
push!(lines, "  e/v entrant share of the vacancy stock (= entrant_vac_share, steady_state.jl:428);")
push!(lines, "  v_pre/v predetermined (inherited) share, exposed to the destruction shock;")
push!(lines, "  T threshold: rise in log e per UNIT rise in the monthly delta_e (level), i.e. e-hat/d.")
push!(lines, "    Per 1 percent rise in delta_e the threshold is delta_e*T = v_pre/(v+L). See beveridge_free_entry.jl.")
push!(lines, "")
push!(lines, "IDENTITIES (verified by assertion above):")
push!(lines, "  e/v     = (delta_e/tau) / kappa_v")
push!(lines, "  v_pre/v = 1 - (delta_e/tau)/kappa_v")
push!(lines, "  T       = kappa_v/(delta_e/tau) - 1")
push!(lines, "")
for c in CALIBS
    b = ss_block(c.delta_e)
    push!(lines, @sprintf("  %-40s s = %.6f   %s", c.name, b.s, c.note))
end
push!(lines, "")
push!(lines, "READING")
push!(lines, "-"^60)
bD1 = ss_block(delta_e_from_ann(0.0320))
bCK = ss_block(SEP)
push!(lines, @sprintf("The threshold is %.1f times higher at the D1 calibration than at Coles-Kelishomi's",
                      (bD1.v_pre / bD1.e) / (bCK.v_pre / bCK.e)))
push!(lines, @sprintf("(%.2f vs %.3f), because delta_e = tau makes the standing vacancy stock %.0f%% new",
                      bD1.v_pre / bD1.e, bCK.v_pre / bCK.e, 100 * bCK.e / bCK.v))
push!(lines, @sprintf("postings, leaving only %.0f%% inherited for a destruction shock to remove. At the D1",
                      100 * bCK.v_pre / bCK.v))
push!(lines, @sprintf("value %.0f%% of the stock is inherited.", 100 * bD1.v_pre / bD1.v))
push!(lines, "")
push!(lines, "delta_e = tau is the maximum admissible destruction rate (s >= 0), so CK sit at the")
push!(lines, "point that MINIMIZES the entry threshold over the whole admissible range.")
push!(lines, "="^104)

text = join(lines, "\n")
println(text)
open("predetermined_vacancy_share.txt", "w") do io
    write(io, text * "\n")
end
println("\nWritten: predetermined_vacancy_share.txt")
