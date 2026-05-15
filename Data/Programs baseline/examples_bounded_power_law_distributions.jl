using Plots

μ = 0.75
b = 1.0

alphas = [1.0, 0.75, 1.25]
labels = ["α=1 (Uniform)", "α=0.75", "α=1.25"]
colors = [:blue, :orange, :green]

# Lower limit is rescaled
a1 = μ*(alphas[1]+1) - alphas[1]
a2 = μ*(alphas[2]+1) - alphas[2]
a3 = μ*(alphas[3]+1) - alphas[3]
a_vals = [a1, a2, a3]

# PDF function generator
function powerlaw_pdf(α, a, b)
    norm = α / (b - a)^α
    return x -> (x >= a && x <= b) ? norm * (x - a)^(α-1) : 0.0
end

pdfs = [powerlaw_pdf(alphas[i], a_vals[i], b) for i in 1:3]
xs = range(minimum(a_vals), 1, length=500)
ys = [[pdf(x) for x in xs] for pdf in pdfs]

plot()
for i in 1:3
    plot!(xs, ys[i], label=labels[i]*", a=$(round(a_vals[i], digits=3))", lw=2, color=colors[i])
end
xlabel!("x")
ylabel!("density")
title!("Bounded Power Law Distributions with Mean 0.75")
