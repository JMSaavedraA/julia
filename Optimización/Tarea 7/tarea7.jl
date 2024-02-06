# Jose Miguel Saavedra Aguilar
# Master in Applied Mathematics
# CIMAT. Optimization. Homework 7

include("trainMNIST.jl")
using Random, Plots, LaTeXStrings

Random.seed!(2023)
β0 = 2*rand(785) .- 1;

# Fit the model
βPF, gPF, kPF, FPF, GPF = pasoFijo(β0,2000,1e-8,1e-1,500); # Stochastic Gradient Method
βN, gN, kN, FN, GN = descensoNewton(β0,1000,1e-8,1e-2,500); # Stochastic modified Newton's Method

# Plot the gradients and save the plot as PNG
plt1 = plot(1:2001,GPF,label = L"\|\|\nabla h(z_t,x^{(i_k)}) \|\|",title="Stochastic Gradient Method");
savefig(plt1,"gradientSGM.png");
plt1 = plot(1:1001,GN,label = L"\|\|\nabla h(z_t,x^{(i_k)}) \|\|",title="Stochastic Modified Newton's Method");
savefig(plt1,"gradientSMNM.png");


# Predict for the test subsample
guessPF = guess(βPF,test_x);
guessN = guess(βN,test_x);

# Calculate the prediction error
errorPF = guessError(guessPF,test_ones);
errorN = guessError(guessN,test_ones);