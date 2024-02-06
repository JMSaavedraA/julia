#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics
#Optimization. Homework 8

using LinearAlgebra
include("trustRegion.jl")
include("funcionesPrueba.jl")
include("filterTrustRegion.jl")
include("optimizadores.jl")


function testTarea(f::Function,∇::Function,H::Function,xOpt::AbstractVector,iter::Int)
    # Execute iter instances of Trust-Region with naive and dogleg from random starting points
    # Initialize:
    avN = 0.0;
    avD = 0.0;
    avB = 0.0;
    avF = 0.0;
    tN = 0.0;
    tD = 0.0;
    tB = 0.0;
    tF = 0.0;
    n = length(xOpt)
    for i = 1:iter
        # Solve iter problems
        x0 = xOpt + 4*rand(n) .- 2 # Take a random starting point
        tN += @elapsed _, _, kN = trustRegionNaive(f,x0,∇,H);
        tD += @elapsed _, _, kD = trustRegionDogleg(f,x0,∇,H);
        tB += @elapsed _, _, kB = descensoNewton(f,x0,∇,H);
        tF += @elapsed _, _, kF = trustRegionRFTR(f,x0,∇,H);
        avN += kN / iter;
        avD += kD / iter;
        avB += kB / iter;
        avF += kF / iter;
    end
    tN /= iter;
    tD /= iter;
    tB /= iter;
    tF /= iter;
    return avN, avD, avB, avF, tN, tD, tB, tF
end

