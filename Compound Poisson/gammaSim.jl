using Random, Distributions, Statistics, Base.Threads, Plots, LaTeXStrings, SpecialFunctions, Printf

Random.seed!(18092023)
# Define the algorithm parameters
n = 1000   # The number of samples of the compound Poisson process
m = 1000
M = 1000
# Define the algorithm parameters
λ = 1
α = 1e-2     # For α ∈ (1,2), the Pareto distribution has infinite variance and finite mean
θ = 1/α     # So the mean is 1
p = zeros(M)
β = 1

g = gamma(α)
ei = expint(-α,α*β)
μ = (λ+1)*((β * (1 + (α*β)^α * ei/g)) - 1)

# Define the random variables distribution
P = Poisson(λ)
Z = Gamma(α,θ)

@threads for k = 1:M
    z = zeros(m)
    for j = 1:m
        s = zeros(n)
        for i = 1:n
            N = 1 + rand(P)
            x = max.((β .- rand(Z,N)), zeros(N))
            s[i] = sum(x)
        end
        z[j] = count(s .> μ)/n
    end
    p[k] = maximum(z)
end

c = count(p .> (1-exp(-1)))
ν = mean(p)
σ = sqrt(var(p))
Φ = Normal(ν,σ)
Plot1 = histogram(p,normalize=:pdf, color=:red,xlims=(0,1),label = L"\mathbb{P}(A)",title = @sprintf("La media es %0.4f, exceden %i",ν,c),size=(1200,800))
y = LinRange(0,1,500)
plot!(y,pdf.(Φ,y),color=:blue,label=L"\mathcal{N}")

savefig(Plot1,"figura1.png")