# R8_snapshot_pre.jl
# =============================================================================
# Q1 step 1 (Notes/Q1_R8_D1_PATHB_plan.md): snapshot of the PRE-R8 code's impulse
# responses, the regression reference for the timing restructure.
#
# Two calibrations, both at the D1 destruction rate (dest_ann = 0.0320) under
# PATH B settings (b_ratio = 0.9, x_v = 0.5):
#   :p0zero — p_0 = 0, no endogenous exit. Here Λ ≡ 1, so the code's exit dating
#             and the draft's coincide: the restructured code must reproduce these
#             IRFs to machine precision (Test 1).
#   :base   — p_0 = 0.5, ω_δ = 0.5, dest_elast_target = 5. Reference only: the
#             restructure is SUPPOSED to change these (exit dated t+1, not t).
# Shocks: z, δ, s, one standard deviation each, 60 months. The s process uses the
# part6b calibration (ρ_s = 0.874, σ_s = 0.0854), not the 0.90/0.010 placeholder.
#
# Output: R8_snapshot_pre.jls  (Dict of case => (; ss, cal, irf_z, irf_δ, irf_s))
#         R8_snapshot_pre.txt  (steady states and selected IRF values)
# =============================================================================

include("run_solution_core.jl")
using Printf

# Numeric shock processes. Named apart from the symbolic ρ_z, σ_z, ... declared in
# run_solution_core.jl, which gen_model_equations() reads at call time.
const SHK = (ρz = 0.902, σz = 0.0092, ρδ = 0.592, σδ = 0.0669, ρs = 0.874, σs = 0.0854)
const T_IR = 60

base_targets = merge(TARGETS, (dest_ann = 0.0320, b_ratio = 0.9, x_v = 0.5))

cases = Dict(
    :p0zero => merge(base_targets, (p_0 = 0.0, dest_end_frac = 0.0, Xc_Y = 0.0)),
    :base   => merge(base_targets, (p_0 = 0.5, dest_end_frac = 0.5, dest_elast_target = 5.0)),
)

function solve_case(tg)
    cal = calibrate_shares(tg)
    @unpack f_e, δ, s, z, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, x_m, ψ, f_m, p_0 = cal
    PAR = [f_e; z; δ; s; b; ϕ; r; σ; ε; A; η_L; κ; ξ_inv; x_m; ψ; f_m; p_0;
           SHK.ρz; SHK.σz; SHK.ρδ; SHK.σδ; SHK.ρs; SHK.σs]
    model = (parameters = parameters, estimate = estimate, estimation = position,
             npar = length(parameters), ns = length(estimate), priors = priors,
             x = x, y = y, xp = xp, yp = yp, variables = variables, varnames = varnames,
             nx = nx, ny = ny, nvar = nvar, e = ex, eta = eta, ne = ne,
             f = gen_model_equations(), nf = nvar,
             SS = SS_symbolics(parameters, tg), PAR_SS = parameters[:],
             flag_order = flag_order, flag_deviation = flag_deviation,
             flag_SSsolver = flag_SSsolver)
    process_model(model)
    sol = solution_interface(model, PAR, SS_numeric(PAR, tg))

    # One shock at a time: column j of the state impact matrix, by state name
    m1 = (; model..., ne = 1)
    function irf(state, sd)
        col = zeros(nx, 1); col[findfirst(==(state), Symbol.(x)), 1] = sd
        100 .* DataFrame(simulate_model(m1, sol.sol_mat, T_IR, col, sol.SS, true, true), varnames)
    end
    return (; ss = sol.ss, cal, irf_z = irf(:z, SHK.σz), irf_δ = irf(:δ, SHK.σδ), irf_s = irf(:s, SHK.σs))
end

results = Dict(k => solve_case(tg) for (k, tg) in cases)
serialize("R8_snapshot_pre.jls", results)

open("R8_snapshot_pre.txt", "w") do io
    println(io, "R8 pre-restructure snapshot — ", length(results), " cases, T = $T_IR months, 100 × log deviations")
    for (k, r) in sort(collect(results); by = first)
        println(io, "\n== case $k ==")
        @printf(io, "  SS: u = %.6f  v = %.6f  θ = %.6f  N = %.6f  δ_e = %.6f  x_c = %.6f\n",
                r.ss.u, r.ss.v, r.ss.θ, r.ss.N, r.ss.δ_e, r.ss.x_c)
        for (lab, df) in (("z", r.irf_z), ("δ", r.irf_δ), ("s", r.irf_s))
            # simulate_model enters the shock in the states at row 2 (row 1 is zero)
            @printf(io, "  %s shock: u[2,3,13] = %.6f %.6f %.6f   v[2,3,13] = %.6f %.6f %.6f   N[13] = %.6f\n",
                    lab, df.u[2], df.u[3], df.u[13], df.v[2], df.v[3], df.v[13], df.N[13])
        end
    end
end
println(read("R8_snapshot_pre.txt", String))
