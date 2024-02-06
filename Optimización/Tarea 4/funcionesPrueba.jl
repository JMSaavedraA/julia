#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics
#Optimization

using Random, DelimitedFiles, LinearAlgebra

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

function Rosenbrock(x::AbstractVector)
    #Rosenbrock function for any size
    n = length(x)
    y = sum(100*(x[2:n] .- x[1:n-1].^2).^2 + (1 .- x[1:n-1]).^2)
    return y
end

function gradienteRosenbrock(x::AbstractVector)
    #Rosenbrock function's gradient for any size
    n = length(x)
    g = zeros(Float64,n)
    g[1] = 400*x[1]*(x[1]^2 - x[2]) + 2*(x[1]-1)
    g[2:n-1] = 200*(x[2:n-1] .- x[1:n-2].^2) + 400*x[2:n-1].*(x[2:n-1].^2 .- x[3:n]) + 2*(x[2:n-1] .- 1)
    g[n] = 200*(x[n] - x[n-1]^2)
    return g
end

function hessianaRosenbrock(x::AbstractVector)
    #Rosenbrock function's hessian matrix for any size
    n = length(x)
    dv = zeros(Float64,n)
    ev = zeros(Float64,n-1)
    dv[1] = -400*(x[2] - 3*x[1]^2) + 2
    dv[2:n-1] = 200 .- 400*(x[3:n] .- 3*x[2:n-1].^2) .+ 2
    dv[n] = 200
    ev .= -400*x[1:n-1]
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

function funImagen(x::AbstractVector)
    #The function in part (3) of HW 4
    n = length(x)
    y = norm(x-Y)^2 + λ*norm(x[2:n]-x[1:n-1])^2
    return y
end

function gradImagen(x::AbstractVector)
    #The function in part (3) of HW 4's gradient
    n = length(x)
    g = zeros(n)
    g[1] = 2*(1+λ)*x[1] - 2*(Y[1] + λ*x[2])
    g[2:n-1] = 2*(1+2*λ)*x[2:n-1] - 2*(Y[2:n-1] + λ*(x[3:n]+x[1:n-2]))
    g[n] = 2*(1+λ)*x[n] - 2*(Y[n] + λ*x[n-1])
    return g
end

function hessianaImagen(x::AbstractVector)
    #The function in part (3) of HW 4's hessian matrix
    n = length(x)
    dv = 2*(1+2*λ)*ones(n)
    ev = -2*λ*ones(n-1)
    dv[1] = 2*(1+λ)
    dv[n] = 2*(1+λ)
    H = SymTridiagonal(dv, ev) #We only store the two diagonals for the whole matrix. This is way faster to solve the system each iterate
    return H
end

