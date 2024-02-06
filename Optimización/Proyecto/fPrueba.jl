#Jose Miguel Saavedra Aguilar and Carlos Manuel Bosch Machado
#CIMAT Master in Applied Mathematics
#Optimization. Final Project

include("trustRegion.jl")
include("funcionesPrueba.jl")
Random.seed!(2023)

x0 = zeros(100);
xR, NgR, kR, ϵR = trustRegionDogleg(Rosenbrock100,x0)
x0 = zeros(4);
xW, NgW, kW, ϵW = trustRegionDogleg(Wood,x0)
x0 = zeros(2);
xR2, NgR2, kR2, ϵR2 = trustRegionDogleg(Rosenbrock2,x0)