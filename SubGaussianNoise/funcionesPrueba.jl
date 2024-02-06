#Jose Miguel Saavedra Aguilar and Carlos Manuel Bosch Machado
#CIMAT Master in Applied Mathematics
#Optimization. Final Project

using LinearAlgebra, Random

function Rosenbrock2(x::AbstractVector)
    #Rosenbrock function for n=2 with random Normal Noise
    y = 100*(x[2]-x[1]^2)^2 + (1-x[1])^2 + 1e-4*randn();
    return y
end

function Rosenbrock100(x::AbstractVector)
    #Rosenbrock function for n=100 with random Normal Noise
    y = sum(100*(x[2:100] .- x[1:99].^2).^2 + (1 .- x[1:99]).^2) + 1e-6*randn();
    return y
end

function Wood(x::AbstractVector)
    #Wood (Colville) function with random Uniform Noise
    y = 100*(x[1]^2 - x[2])^2 + (x[1] - 1)^2 + (x[3] - 1)^2 + 90*(x[3]^2 - x[4])^2 + 10.1*((x[2] - 1)^2 + (x[4] - 1)^2) + 19.8*(x[2] - 1)*(x[4] - 1) + 1e-5*rand();
    return y
end
