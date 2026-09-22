# run_xi_eps_sweep.jl
# Joint sweep: ξ_inv × ε to find if ANY combination delays δ→u peak.
# BED spec (dest_ann=0.0320), everything else at defaults.

include("run_solution_core.jl")

ρ_z = 0.902;  σ_z = 0.0092
ρ_δ = 0.592;  σ_δ = 0.0669
ρ_s = 0.8741; σ_s = 0.0854

T_IR    = 63
flag_IR = true; flag_logdev = true

include("time_series_fun.jl")

function solve_one(model, xi_val, eps_val)
    tgt = (TARGETS...,
        dest_ann          = 0.0320,
        dest_end_frac     = 0.5,
        p_0               = 0.5,
        dest_elast_target = 5.0,
        b_ratio           = 0.9,
        x_v               = 0.5,
        ξ_inv             = xi_val,
        ε                 = eps_val,
    )
    cal = calibrate_shares(tgt)
    @unpack f_e, δ, s, z, b, ϕ, r, σ, ε, A, η_L, κ, ξ_inv, x_m, ψ, f_m, p_0 = cal
    s = max(s, 0.0)
    PAR = [f_e; z; δ; s; b; ϕ; r; σ; ε; A; η_L; κ; ξ_inv; x_m; ψ; f_m; p_0;
           ρ_z; σ_z; ρ_δ; σ_δ; ρ_s; σ_s]
    SS_val = SS_numeric(PAR, tgt)
    sol = solution_interface(model, PAR, SS_val)
    @unpack ss, SS, sol_mat = sol

    eta_δ_col = reshape([0.0, 0.0, 0.0, 0.0, σ_δ, 0.0], 6, 1)
    model_1 = (; model..., ne = 1)
    irf_δ = 100 .* DataFrame(
        simulate_model(model_1, sol_mat, T_IR, eta_δ_col, SS, flag_IR, flag_logdev),
        varnames)

    to_pp(col, ss_val) = (exp.(col ./ 100) .- 1) .* ss_val .* 100
    u_pp = to_pp(irf_δ.u, ss.u)
    q_avg(x) = [mean(x[i:i+2]) for i in 1:3:length(x)-2]
    u_q = q_avg(u_pp)
    peak_h = argmax(abs.(u_q)) - 1

    # Also report the hx eigenvalues to understand dynamics
    hx = sol.sol_mat.hx
    eigvals_hx = abs.(eigvals(hx))
    max_eig = maximum(eigvals_hx)

    return (xi=xi_val, eps=eps_val, u_peak_h=peak_h,
            u_peak=round(u_q[peak_h+1], sigdigits=3),
            u_h20=round(u_q[end], sigdigits=3),
            max_eig=round(max_eig, sigdigits=4),
            u_q=u_q)
end

using LinearAlgebra

println("Building model...")
dummy_targets = (TARGETS..., dest_ann=0.0320, dest_end_frac=0.5, p_0=0.5,
                 dest_elast_target=5.0, b_ratio=0.9, x_v=0.5, ξ_inv=1.0)
SS_sym = SS_symbolics(parameters, dummy_targets)
model = (parameters=parameters, estimate=estimate, estimation=position,
         npar=length(parameters), ns=length(estimate), priors=priors,
         x=x, y=y, xp=xp, yp=yp, variables=variables, varnames=varnames,
         nx=nx, ny=ny, nvar=nvar, e=ex, eta=eta, ne=ne,
         f=gen_model_equations(), nf=nvar,
         SS=SS_sym, PAR_SS=parameters[:],
         flag_order=flag_order, flag_deviation=flag_deviation,
         flag_SSsolver=flag_SSsolver)
process_model(model)
println("Model ready.\n")

xi_vals  = [1.0, 3.0, 5.0, 8.0]
eps_vals = [4.3, 3.0, 2.0, 1.5]

println("ε        ξ_inv    peak_h   u_peak(pp)  u_h20(pp)   ratio    max|eig(hx)|")
println("-"^80)
for eps in eps_vals
    for xi in xi_vals
        try
            r = solve_one(model, xi, eps)
            ratio = round(r.u_h20 / r.u_peak, sigdigits=3)
            println("$(rpad(eps, 9))$(rpad(xi, 9))$(rpad(r.u_peak_h, 9))$(rpad(r.u_peak, 12))$(rpad(r.u_h20, 12))$(rpad(ratio, 9))$(r.max_eig)")
        catch ex
            println("$(rpad(eps, 9))$(rpad(xi, 9)) FAILED: $(sprint(showerror, ex))")
        end
    end
end

# Also print the full quarterly IRF for the most extreme case
println("\n\nFull quarterly δ→u IRF (pp) for most extreme case (ξ_inv=8, ε=1.5):")
try
    r = solve_one(model, 8.0, 1.5)
    for (h, val) in enumerate(r.u_q)
        println("  h=$(rpad(h-1, 4)) $(round(val, sigdigits=4))")
    end
catch ex
    println("  FAILED: $(sprint(showerror, ex))")
end
