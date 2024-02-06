#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics
#Optimization

using LinearAlgebra

function Rosenbrock2(x::AbstractVector)
    #Rosenbrock function for n=2
    y = 100*(x[2]-x[1]^2)^2 + (1-x[1])^2
    return y
end

function gradienteRosenbrock2(x::AbstractVector)
    #Rosenbrock function's gradient for n=2
    g = zeros(Float64,2)
    g[1] = 400*x[1]*(x[1]^2 - x[2]) + 2*(x[1]-1)
    g[2] = 200*(x[2] - x[1]^2)
    return g
end

function hessianaRosenbrock2(x::AbstractVector)
    #Rosenbrock function's hessian matrix for n=2
    H = zeros(Float64,2,2)
    H[1,1] = -400*(x[2] - 3*x[1]^2) + 2
    H[2,1] = -400*x[1]
    H[1,2] = -400*x[1]
    H[2,2] = 200
    return H
end

