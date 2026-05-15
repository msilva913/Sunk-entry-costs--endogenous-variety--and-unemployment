# --- Synthetic Hill estimator demo in Julia ---
using Plots
using Printf
# Simple Pareto generator (xmin > 0, alpha > 0)
# Inverse-CDF method: X = xmin * (1 - U)^(-1/alpha)
function rand_pareto(n::Integer; xmin::Float64=1.0, alpha::Float64=2.0)
    u = rand(n)
    return xmin * (1 .- u) .^ (-1/alpha)
end

# Optionally build a mixed distribution: many small observations + Pareto tail
function rand_mixed(n::Integer; frac_tail::Float64=0.05, xmin::Float64=1.0, alpha::Float64=2.0)
    n_tail = max(1, round(Int, n * frac_tail))
    n_body = n - n_tail
    # body: lognormal small masses
    body = exp.(randn(n_body) .* 0.6 .+ 0.0)  # adjust for body spread
    tail = rand_pareto(n_tail; xmin=xmin, alpha=alpha)
    return vcat(body, tail)
end

# Hill estimator functions (gamma = 1/alpha)
function hill_gamma(data::AbstractVector{<:Real}, k::Integer)
    x = collect(filter(!isnan, data))
    x = x[x .> 0]            # strictly positive
    n = length(x)
    @assert k >= 1 && k < n "k must be between 1 and n-1"
    #γ_k = (1/k)\sum_{i=1}^k[log X(i) - log X_(k+1)]
    xs = sort(x, rev=true)   # descending order(rev denotes reverse)
    # Calculate #
    threshold = xs[k+1]
    logs = log.(xs[1:k]) .- log(threshold) 
    γ_hat = mean(logs)

    se_γ = γ_hat / sqrt(k)   # asymptotic approx
    return γ_hat, se_γ
end

function hill_alpha(data::AbstractVector{<:Real}, k::Integer)
    γ, seγ = hill_gamma(data, k)
    α_hat = 1.0 / γ
    seα = α_hat / sqrt(k)    # delta method approx: se(α) ≈ α / sqrt(k)
    return α_hat, seα, γ, seγ
end

# Compute Hill estimates across a range of k and return arrays
function hill_series(data::AbstractVector{<:Real}; kmin::Integer=5, kmax::Union{Nothing,Int}=nothing)
    x = collect(filter(!isnan, data))
    x = x[x .> 0]
    n = length(x)

    kmax = isnothing(kmax) ? floor(Int, n/2) : min(kmax, n-1)
    kmin = max(1, min(kmin, kmax-1))
    ks = kmin:kmax
    αs = Float64[]
    seαs = Float64[]
    γs = Float64[]
    for k in ks
        α_hat, seα, γ_hat, seγ = hill_alpha(x, k)
        push!(αs, α_hat)
        push!(seαs, seα)
        push!(γs, γ_hat)
    end
    return ks, αs, seαs, γs
end

# Optional plotting helper (requires Plots.jl)




# ------------------ Demo run ------------------

# 1) Pure Pareto demo
n = 10000
true_alpha = 2.0       # α_true
data_pareto = rand_pareto(n; xmin=1.0, alpha=true_alpha)

# choose some k values (top fractions). Example: top 0.5%, 1%, 2%, 5%
ks = [max(5, floor(Int, p * n)) for p in (0.005, 0.01, 0.02, 0.05)]
println("Pure Pareto (α_true = $(true_alpha)), n = $n")
for k in ks
    α_hat, seα, γ_hat, seγ = hill_alpha(data_pareto, k)
    println(" k=$(k) (top $(round(100*k/n,digits=3))%): α̂=$(round(α_hat,digits=4)) se(α̂)=$(round(seα,digits=4)), γ̂=$(round(γ_hat,digits=4)) se(γ̂)=$(round(seγ,digits=4))")
end

# 2) Mixed body + Pareto tail demo (more realistic)
data_mixed = rand_mixed(n; frac_tail=0.05, xmin=1.0, alpha=true_alpha)
println("\nMixed distribution (5% tail) demo, n = $n")
ks_range = floor.(Int, range(max(5, round(Int,0.002*n)), stop=floor(Int,0.08*n), length=20))
ks_range = unique(sort(ks_range))
ks_range = ks_range[ks_range .< n]  # safe
ks_range = ks_range[1:min(end,50)]

# compute series and print a short table
println(" k   α̂    se(α̂)   γ̂   se(γ̂)")
for k in ks_range
    α_hat, seα, γ_hat, seγ = hill_alpha(data_mixed, k)
    println(@sprintf("%4d %6.3f %8.3f %6.3f %8.3f", k, α_hat, seα, γ_hat, seγ))
end

# Produce a hill plot for the mixed data (γ̂ vs k)

ks_plot, αs, seαs, γs = hill_series(data_mixed; kmin=10, kmax=floor(Int,n*0.08))
# plot γ̂(k)


Plots.plot(ks, vals; xlabel="k (top order stats used)", ylabel="γ̂ (Hill)", title="Hill plot: γ̂ vs k (mixed data)", legend=false)
# If you want α̂ plot instead, run:
# hill_plot(ks_plot, αs; ylabel="α̂", title="Hill plot: α̂ vs k (mixed data)")

# End of demo