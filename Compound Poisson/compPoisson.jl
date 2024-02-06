using Distributions, Statistics, Base.Threads, Plots, LaTeXStrings

# Define the algorithm parameters
n = 1000   # The number of samples of the compound Poisson process
m = 1000
α = 3/2     # For α ∈ (1,2), the Pareto distribution has infinite variance and finite mean
θ = 1/3     # So the mean is 1
λ = 10
p = zeros(m)


# Define the random variables distribution
P = Poisson(λ)
Z = Pareto(α,θ)

@threads for j = 1:m
    s = zeros(n)
    for i = 1:n
        N = 1 + rand(P)
        x = rand(Z,N)
        s[i] = sum(x)
    end
    p[j] = count(s .> (λ+1))/n
end


histogram(p,normalize=:pdf, color=:red,xlims=(0,1),label = L"\mathbb{P}(A)")
μ = mean(p)
σ = sqrt(var(p))
Φ = Normal(μ,σ)
y = LinRange(0,1,500)
plot!(y, pdf.(Φ,y), color=:blue, label = L"\mathcal{N}")