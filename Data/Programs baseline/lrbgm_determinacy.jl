# lrbgm_determinacy.jl
# =============================================================================
# Determinacy of the LR-BGM business formation block (rem:bgm_nest) as ε varies.
# p_0 = 0 (X^c = 0), exogenous constant L and z, DS-CES ρ(N) = N^{1/(ε−1)}, μ = ε/(ε−1).
# One state N_t, one jump: entry labor x_t = L^e_t.
#   N_{t+1} = (1−δ)(N_t + z x_t / f_e)                                      (eq:N_lom_eq)
#   f_e ρ(N_t) = β(1−δ)(C_t/C_{t+1})^σ {f_e ρ(N_{t+1}) + (μ−1) C_{t+1}/N_{t+1}}
#                                       with C_t = ρ(N_t) z (L − x_t)        (eq:euler_bf)
# Dividing by f_e ρ(N_{t+1}): the dividend (μ−1) z L^c/(f_e N) is free of ρ (dilution only), but the
# capital gain ρ(N_{t+1})/ρ(N_t) and C_t/C_{t+1} keep a 1/(ε−1) elasticity. Determinacy
# needs exactly one generalized eigenvalue outside the unit circle. f_e normalizes N̄ = 1.
# Output: lrbgm_determinacy.txt
# =============================================================================
using LinearAlgebra, Printf

function lrbgm_roots(ε, σ; δe = 1 - (1 - 0.0320)^(1/12), r = 1.04^(1/12) - 1, L = 0.93, z = 1.0)
    β, μ, δ = 1 / (1 + r), ε / (ε - 1), δe          # p_0 = 0 ⇒ δ_e = δ
    ρf(N) = N^(1 / (ε - 1))
    # steady state with N̄ = 1: x̄ = δ f_e/((1−δ) z); f_e(1−β(1−δ)) = β(1−δ)(μ−1) z (L − x̄)
    k = β * (1 - δ) * (μ - 1) * z
    fe = k * L / ((1 - β * (1 - δ)) + k * δ / ((1 - δ) * z))
    xbar = δ * fe / ((1 - δ) * z)
    lom(N, x, Np) = Np - (1 - δ) * (N + z * x / fe)
    C(N, x) = ρf(N) * z * (L - x)
    eul(N, x, Np, xp) = fe * ρf(N) - β * (1 - δ) * (C(N, x) / C(Np, xp))^σ *
                        (fe * ρf(Np) + (μ - 1) * C(Np, xp) / Np)
    s = [1.0, xbar]                                  # (N, x) at steady state
    F(v) = [lom(v[1], v[2], v[3]), eul(v[1], v[2], v[3], v[4])]
    @assert maximum(abs, F([s; s])) < 1e-12
    h = 1e-6
    J = hcat([(F([s; s] .+ h .* (1:4 .== j)) .- F([s; s] .- h .* (1:4 .== j))) ./ (2h) for j in 1:4]...)
    A, B = J[:, 3:4], -J[:, 1:2]                     # A [N', x'] = B [N, x]
    λ = eigvals(B, A)
    return sort(abs.(λ))
end

lines = String[]
log!(s) = (push!(lines, s); println(s))
log!("LR-BGM block (p_0 = 0, exogenous L): |generalized eigenvalues|; determinate iff exactly one > 1")
for σ in (1.0, 2.0), ε in (1.05, 1.2, 1.4, 1.5, 2.0, 2.5, 3.0, 4.3, 8.0)
    λ = lrbgm_roots(ε, σ)
    nout = count(>(1), λ)
    log!(@sprintf("  σ = %.1f  ε = %4.2f   |λ| = %.6f, %.4g   %s", σ, ε, λ[1], λ[2],
                  nout == 1 ? "determinate" : nout == 0 ? "INDETERMINATE" : "NO bounded equilibrium"))
end
write("lrbgm_determinacy.txt", join(lines, "\n") * "\n")
