# gs_finite_xi_entry_offset.jl
# =============================================================================
# The δ shock at a FINITE entry elasticity ξ in the economy of prop:fe (AGS with
# σ = κ = 0 = GS, p_0 = 0, costless reposting), draft timing (S10). Task R12c: decides
# the form of the finite-ξ sufficient condition in the general-ξ lemma (old Prop. 5
# Part 3). At finite ξ the system is not block recursive: u and v_pre are states.
#
# Method. The linearized perfect-foresight path to a known AR(1) shock path equals the
# first-order rational-expectations IRF (certainty equivalence). The model is stacked
# over t = 0..T, the Jacobian is taken by central differences around the steady state,
# and the linear system is solved with `\`. Robustness to T is asserted.
#
#   u_{t+1}  = [1 − (1−δ_t) f(θ_t)] u_t + [δ_t + s(1−δ_t)](1 − u_t)
#   vp_{t+1} = (1−δ_t)[(1 − q(θ_t)) v_t + s(1 − u_t)]
#   K_t      = Q_t − β(1−δ_t) Q_{t+1}                                     eq:Kdef
#   κ + K_t/q_t = β(1−δ_t)[(1−ϕ)(y − K_{t+1} − b) − ϕθ_{t+1}(K_{t+1} + q_{t+1}κ)
#                          + (1−s)(κ + K_{t+1}/q_{t+1})]                   eq:jcc_eq
#   v_t = vp_t + e(Q_t),  e(Q) = (Q/x_m)^ξ,  θ_t = v_t/u_t,  δ_t = δ̄ exp(δ̃_t), δ̃_t = ρ^t δ̃_0
#
# Entry offset at h = 1 (old Part 3): destruction of the stock carried into t = 1 is
# X̄·δ̄·δ̃_0 with X̄ = (1−q̄)v̄ + s̄(1−ū); entry adds (1−δ̄)[1−(1−η)q̄]Δe_0 + Δe_1.
#   offset ≡ entry contribution / destruction;  v_1 < v̄  iff  offset < 1.
#
# OUTPUT  gs_finite_xi_entry_offset.txt
# =============================================================================

include(joinpath(@__DIR__, "beveridge_free_entry.jl"))   # calibrate, CALIBS, TG, ρ_δ
using Printf, LinearAlgebra, SparseArrays

"""
Linear IRF to δ̃_0 = 1 (log shock); level deviations per unit δ̃_0. Unknowns u, v_pre,
ℓ = ln e and K. With Q = x_m exp(ℓ/ξ), eq:Kdef is written in deviations through expm1,
so the system stays well conditioned as ξ grows.
"""
function irf_xi(c, ξ; ρ = ρ_δ, T = 2000, nl_check = false)
    η, s, ϕ, b, κ, y = c.η, c.s, c.ϕ, c.b, c.κ, c.w_int
    β = 1 / (1 + c.r);  A = c.q * c.θ^c.η
    q(θ) = A * θ^(-η);  f(θ) = A * θ^(1 - η)
    ℓ̄ = log(c.e);  x_m = c.Q / exp(ℓ̄ / ξ)
    dQ(ℓ) = x_m * (expm1(ℓ / ξ) - expm1(ℓ̄ / ξ))        # Q(ℓ) − Q̄
    vp̄ = c.v - c.e
    n = T + 1
    function resid(X, d0)
        u = X[1:n]; vp = X[n+1:2n]; ℓ = X[2n+1:3n]; K = X[3n+1:4n]
        δ = [c.δ_e * exp(d0 * ρ^t) for t in 0:T]
        v = vp .+ exp.(ℓ);  θ = v ./ u
        R = zeros(4n)
        R[1] = u[1] - c.u;  R[2] = vp[1] - vp̄
        for t in 1:T                                   # laws of motion, t-1 → t (0-based)
            i = t
            R[2+t]   = u[i+1] - ((1 - (1 - δ[i]) * f(θ[i])) * u[i] + (δ[i] + s * (1 - δ[i])) * (1 - u[i]))
            R[2+T+t] = vp[i+1] - (1 - δ[i]) * ((1 - q(θ[i])) * v[i] + s * (1 - u[i]))
        end
        for i in 1:n                                   # forward-looking block
            dQn = i < n ? dQ(ℓ[i+1]) : 0.0
            Kn = i < n ? K[i+1] : c.K
            θn = i < n ? θ[i+1] : c.θ
            # K − K̄ = (Q − Q̄) − β(1−δ)(Q' − Q̄) + β(δ − δ̄)Q̄
            R[2+2T+i] = (K[i] - c.K) - dQ(ℓ[i]) + β * (1 - δ[i]) * dQn - β * (δ[i] - c.δ_e) * c.Q
            S = (1 - ϕ) * (y - Kn - b) - ϕ * θn * (Kn + q(θn) * κ) + (1 - s) * (κ + Kn / q(θn))
            R[2+2T+n+i] = κ + K[i] / q(θ[i]) - β * (1 - δ[i]) * S
        end
        return R
    end
    X̄ = vcat(fill(c.u, n), fill(vp̄, n), fill(ℓ̄, n), fill(c.K, n))
    @assert maximum(abs, resid(X̄, 0.0)) < 1e-10 "steady state fails"
    # Sparse Jacobian by period coloring: a residual dated p loads on periods p−1, p, p+1,
    # so perturbing every third period of one block at once identifies each column.
    per = vcat(0, 0, 1:T, 1:T, 0:T, 0:T)               # 0-based period of each residual
    I_, J_, V_ = Int[], Int[], Float64[]
    for blk in 0:3, k in 0:2
        h = 1e-6 * max(abs(X̄[blk*n+1]), 1e-3)
        cols = [blk * n + i + 1 for i in 0:T if i % 3 == k]
        Xp = copy(X̄); Xp[cols] .+= h;  Xm = copy(X̄); Xm[cols] .-= h
        d = (resid(Xp, 0.0) .- resid(Xm, 0.0)) ./ (2h)
        for r in eachindex(d)
            d[r] == 0 && continue
            i = findfirst(j -> 0 <= j <= T && j % 3 == k, per[r]-1:per[r]+1) + per[r] - 2
            push!(I_, r); push!(J_, blk * n + i + 1); push!(V_, d[r])
        end
    end
    Js = sparse(I_, J_, V_, 4n, 4n)
    # Blanchard–Kahn at an interior date p: rows LoM(p+1) and forward(p) give
    # A x_{p+1} = B x_p with x = (u, v_pre, ℓ, K); u and v_pre are the only states.
    p = T ÷ 2;  rows = [2 + p + 1, 2 + T + p + 1, 2 + 2T + p + 1, 2 + 2T + n + p + 1]
    cp = [b * n + p + 1 for b in 0:3]
    Am = Matrix(Js[rows, cp .+ 1]);  Bm = -Matrix(Js[rows, cp])
    nstable = count(λ -> abs(λ) < 1 - 1e-10, eigvals(Bm, Am))
    @assert nstable == 2 "BK fails: $nstable stable roots, ξ = $ξ, δ̄ = $(c.δ_e)"
    J = lu(Js)
    hd = 1e-6
    Rd = (resid(X̄, hd) .- resid(X̄, -hd)) ./ (2hd)
    dX = -(J \ Rd)
    if nl_check                                        # quasi-Newton on the nonlinear path
        d0 = 1e-5;  X = X̄ .+ d0 .* dX
        for _ in 1:50
            Δ = J \ resid(X, d0);  X .-= Δ
            maximum(abs, Δ) < 1e-15 && break
        end
        @assert maximum(abs, resid(X, d0)) < 1e-12 "nonlinear path does not converge, ξ = $ξ"
        @assert isapprox((X .- X̄) ./ d0, dX; rtol = 1e-3, atol = 1e-3 * maximum(abs, dX)) "linear ≠ nonlinear, ξ = $ξ"
    end
    du = dX[1:n]; dvp = dX[n+1:2n]; de = c.e .* dX[2n+1:3n]     # Δe = ē·Δℓ
    dv = dvp .+ de
    dθ = (dv .- c.θ .* du) ./ c.u
    X̄d = (1 - c.q) * c.v + c.s * (1 - c.u)
    destr = X̄d * c.δ_e
    entry = (1 - c.δ_e) * (1 - (1 - η) * c.q) * de[1] + de[2]
    @assert abs(dv[2] + destr - entry) <= 1e-5 * (destr + abs(entry)) "h = 1 decomposition fails, ξ = $ξ, δ̄ = $(c.δ_e)"
    return (; nstable, du, dv, de, dθ, θh = dθ ./ c.θ, offset = entry / destr, v1 = dv[2] / c.v,
              u1 = du[2] / c.u)
end

lines = String[]
P(s) = push!(lines, s)
fmtrow(f, a...) = Printf.format(Printf.Format(f), a...)
P("="^100)
P("δ SHOCK AT FINITE ξ (GS economy of prop:fe, κ = 0): ENTRY OFFSET OF THE h = 1 DESTRUCTION")
P("offset = entry contribution to v_1 / destruction of the stock carried into t = 1;  v̂_1 < 0 iff offset < 1.")
P("="^100)

tg = (TG..., x_v = 1.0)
c1 = calibrate(CALIBS[1][2]; t = tg)

# Check 1: as ξ grows the finite-ξ tightness response converges to prop:fe (log units)
for ρ in (0.0, ρ_δ)
    aF = c1.δ_e * tightness(c1; ρ = ρ).a              # prop:fe c_δ, per unit log shock
    gap = [abs(irf_xi(c1, x; ρ = ρ).θh[1] / aF - 1) for x in (1e4, 1e6, 1e8)]
    @assert gap[1] > gap[2] > gap[3] && gap[3] < 1e-3 "ξ → ∞ limit fails at ρ = $ρ: $gap"
    P(@sprintf("  ξ → ∞ check, ρ = %.3f: |θ̂_0/c_δ − 1| = %.3f, %.3f, %.5f at ξ = 1e4, 1e6, 1e8",
               ρ, gap...))
end
# Check 2: horizon robustness
@assert isapprox(irf_xi(c1, 1.0; T = 1000).offset, irf_xi(c1, 1.0; T = 2000).offset; rtol = 1e-6)

ξs = (0.265, 0.5, 1.0, 2.0, 10.0, 100.0, 1e4)
for ρ in (0.0, ρ_δ)
    P("")
    P(@sprintf("── ρ = %.3f: offset at h = 1 (and v̂_1 per 1%% δ shock, %%) ", ρ) * "─"^56)
    P(fmtrow("%-26s %8s" * " %10s"^length(ξs), "Calibration", "δ̄/τ̄", [@sprintf("ξ=%g", x) for x in ξs]...))
    P("-"^100)
    for (nm, d) in CALIBS
        c = calibrate(d; t = tg)
        Rs = [irf_xi(c, x; ρ = ρ) for x in ξs]
        P(fmtrow("%-26s %8.3f" * " %10.3f"^length(ξs), nm, c.δ_e / c.τ, [R.offset for R in Rs]...))
        P(fmtrow("%-26s %8s" * " %10.3f"^length(ξs), "   v̂_1 (%)", "", [R.v1 for R in Rs]...))
    end
    P("-"^100)
end

# Scaling in δ̄. calibrate() re-targets f̄, q̄, τ̄ and X/Y at every δ̄, so this is a
# cross-economy comparison, not a ceteris paribus limit: as δ̄ → 0 the vacancy asset Q̄
# grows like 1/(r+δ̄) and ϕ falls (printed), and the offset need not vanish.
P("")
P("── Offset across δ̄ (recalibrated economies, κ = 0) ───────────────────────────────────────────────")
P(@sprintf("%-8s %7s %7s %7s %6s | %9s %9s %9s %9s | %9s %9s", "δ̄ ann", "δ̄/τ̄", "ē/v̄", "Q̄", "ϕ",
           "ρ̂ ξ=.265", "ρ̂ ξ=0.5", "ρ̂ ξ=1", "ρ̂ ξ=2", "iid ξ=.265", "iid ξ=1"))
P("-"^100)
a_ok = Float64[]
for a in (0.001, 0.0025, 0.005, 0.01, 0.02, 0.032, 0.05, 0.0754, 0.10, 0.2)
    d = δ_e_from_ann(a)
    d >= TG.sep && continue
    c = try calibrate(d; t = tg) catch; P(@sprintf("%-8.4f  calibration fails", a)); continue end
    o = [irf_xi(c, x).offset for x in (0.265, 0.5, 1.0, 2.0)]
    oi = [irf_xi(c, x; ρ = 0.0).offset for x in (0.265, 1.0)]
    push!(a_ok, a)
    P(@sprintf("%-8.4f %7.4f %7.4f %7.1f %6.3f | %9.3f %9.3f %9.3f %9.3f | %9.3f %9.3f", a,
               c.δ_e / c.τ, c.e / c.v, c.Q, c.ϕ, o..., oi...))
end
P("-"^100)
const a_lo = 0.01         # below ≈1%/yr the recalibrated Q̄ explodes (table above)

# Threshold δ̄*(ξ): the largest δ̄ (annual) with offset < 1, by bisection, ρ = 0.592
function δstar(ξ; ρ = ρ_δ)
    g(a) = irf_xi(calibrate(δ_e_from_ann(a); t = tg), ξ; ρ = ρ).offset - 1
    lo, hi = a_lo, 1 - (1 - TG.sep * 0.999)^12
    sign(g(lo)) == sign(g(hi)) && return (g(hi) < 0 ? Inf : NaN)
    for _ in 1:50
        m = (lo + hi) / 2
        sign(g(m)) == sign(g(lo)) ? (lo = m) : (hi = m)
    end
    return (lo + hi) / 2
end
P("")
P("── δ̄*(ξ): v̂_1 < 0 iff δ̄ < δ̄*(ξ), searched on δ̄ ∈ [1%/yr, τ̄] (annual; all = whole range) ─────────")
P(@sprintf("%-10s %14s %14s", "ξ", "ρ = 0", "ρ = 0.592"))
P("-"^100)
for x in (0.265, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0)
    fmt(z) = isinf(z) ? "all" : isnan(z) ? "none" : @sprintf("%.4f", z)
    P(@sprintf("%-10g %14s %14s", x, fmt(δstar(x; ρ = 0.0)), fmt(δstar(x))))
end
P("-"^100)
P("="^100)

text = join(lines, "\n")
println(text)
open(joinpath(@__DIR__, "gs_finite_xi_entry_offset.txt"), "w") do io
    write(io, text * "\n")
end
