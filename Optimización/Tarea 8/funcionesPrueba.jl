#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics
#Optimization

using LinearAlgebra

function Rosenbrock2(x::AbstractVector)
    #Rosenbrock function for n=2
    y = 100*(x[2]-x[1]^2)^2 + (1-x[1])^2
    return y
end

function Rosenbrock100(x::AbstractVector)
    #Rosenbrock function for n=100
    y = sum(100*(x[2:100] .- x[1:99].^2).^2 + (1 .- x[1:99]).^2)
    return y
end

function gradienteRosenbrock100(x::AbstractVector)
    #Rosenbrock function's gradient for n=100
    g = zeros(Float64,100)
    g[1] = 400*x[1]*(x[1]^2 - x[2]) + 2*(x[1]-1)
    g[2:99] = 200*(x[2:99] .- x[1:98].^2) + 400*x[2:99].*(x[2:99].^2 .- x[3:100]) + 2*(x[2:99] .- 1)
    g[100] = 200*(x[100] - x[99]^2)
    return g
end

function hessianaRosenbrock100(x::AbstractVector)
    #Rosenbrock function's hessian matrix for n=100
    dv = zeros(Float64,100)
    ev = zeros(Float64,99)
    dv[1] = -400*(x[2] - 3*x[1]^2) + 2
    dv[2:99] = 200 .- 400*(x[3:100] .- 3*x[2:99].^2) .+ 2
    dv[100] = 200
    ev .= -400*x[1:99]
    H = SymTridiagonal(dv, ev) #We only store the two diagonals for the whole matrix. This is way faster to solve the system each iterate
    return H
end

function Wood(x::AbstractVector)
    #Wood (Colville) function
    y = 100*(x[1]^2 - x[2])^2 + (x[1] - 1)^2 + (x[3] - 1)^2 + 90*(x[3]^2 - x[4])^2 + 10.1*((x[2] - 1)^2 + (x[4] - 1)^2) + 19.8*(x[2] - 1)*(x[4] - 1)
    return y
end

function gradienteWood(x::AbstractVector)
    #Wood (Colville) function's gradient
    g = zeros(Float64,4)
    g[1] = 400*(x[1]^2 - x[2])*x[1] + 2*(x[1] - 1)
    g[2] = 200*(x[2] - x[1]^2) + 20.2*(x[2] - 1) + 19.8*(x[4] - 1)
    g[3] = 360*(x[3]^2 - x[4])*x[3] + 2*(x[3] - 1)
    g[4] = 180*(x[4] - x[3]^2) + 20.2*(x[4] - 1) + 19.8*(x[2] - 1)
    return g
end

function hessianaWood(x::AbstractVector)
    #Wood (Colville) function's hessian matrix
    H=zeros(Float64,4,4)
    H[1,1] = 400*(3*x[1]^2-x[2]) + 2
    H[1,2] = -400*x[1]
    H[2,2] = 220.2
    H[3,3] = 360*(3*x[3]^2-x[4]) + 2
    H[2,4] = 19.8
    H[3,4] = -360*x[3]
    H[4,4] = 200.2
    sH = Symmetric(H)#We only store half of the matrix this way
    return sH
end

function Branin(x::AbstractVector)
    x1 = x[1]/(2π);
    x2 = x[2];
    f = (x2 +(10 - 5.1*x1 )*x1 -6) + 5*(8π-1)*cospi(2*x1)/(4π) + 10;
    return f
end

function gradienteBranin(x::AbstractVector)
    x1 = x[1]/(2π);
    x2 = x[2];
    g = zeros(2);
    a = 2*(x2 +(10 - 5.1*x1 )*x1 -6);
    g[1] = (a*(5 - 5.1*x1) + 1.25*(1-8π)*sinpi(2*x1))/π;
    g[2] = a
    return g
end

function hessianaBranin(x::AbstractVector)
    x1 = x[1]/(2π);
    x2 = x[2];
    H = zeros(2,2);
    a = x2 +(10 - 5.1*x1 )*x1 -6;
    b = 5 - 5.1*x1;
    H[1,1] = ((2*b^2 - 5.1*a)/π + 1.25*(1-8π)*cospi(2*x1))/π;
    H[1,2] = 2*b/π;
    H[2,1] = 2*b/π;
    H[2,2] = 2;
    H = Symmetric(H)
    return H
end