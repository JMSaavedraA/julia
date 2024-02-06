using Random, Statistics, Distributions, StatsBase, Plots, LaTeXStrings

Random.seed!(2023)
B = 1000
n = 50
k = 20

x = rand(Uniform(),n)
b = zeros(B)

for j = 1:B
    y = sample(x, n)
    b[j] = maximum(y)
end

plt1 = histogram(b, normalize = :pdf, label=L"\hat{\theta}", size=(1600,1200), legend=:topright, palette=:Oranges_3)
f(y) = n*y^(n-1)
t=LinRange(0,1,B)
plot!(t,f, label="Densidad Uniforme", color=:red)
savefig(plt1, "8-7.png")