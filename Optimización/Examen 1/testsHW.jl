#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics
#Optimization.

using LinearAlgebra
include("trustRegion.jl")
include("funcionesPrueba.jl")


function testTarea(f::Function,∇::Function,H::Function,xOpt::AbstractVector,iter::Int)
    # Execute iter instances of Trust-Region with naive and dogleg from random starting points
    # Initialize:
    avD = 0.0;
    avE = 0.0;
    tD = 0.0;
    tE = 0.0;
    CE = 0;
    CD = 0;
    n = length(xOpt)
    for i = 1:iter
        # Solve iter problems
        x0 = xOpt + 4*rand(n) .- 2 # Take a random starting point
        tD += @elapsed _, _, kD, cD = trustRegionDogleg(f,x0,∇,H);
        tE += @elapsed _, _, kE, cE = trustRegionExacta(f,x0,∇,H);
        avD += kD / iter;
        avE += kE / iter;
        CE += cE
        CD += cD
    end
    tE /= iter;
    tD /= iter;
    return avD, avE, tD, tE, CD, CE
end

