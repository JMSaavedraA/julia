using Random
using Printf
include("Integracion.jl")

Random.seed!(2022)
#Definimos la función f_2(x)
function f2(x::AbstractVector)
    z=abs(x[1]+x[2])
    return z
end

#Límites de integración
a=[-1.0, -1.0]
b=[1.0, 1.0]
#valores superiores e inferiores de f_1(x)
fMax=2.0
fMin=0.0
#Calculamos la integral de f_1(x) en D
I=integralMonteCarlo(f2,a,b,10000,1000,1e-15,2,fMax,fMin)
#Error absoluto y relativo
eAbs=abs(I-(8/3))
eRel=eAbs*3/8

println(@sprintf("El valor aproximado de la integral es de %0.7g",I))
println(@sprintf("El error absoluto de integración por el método de Monte Carlo es de %0.7g",eAbs))
println(@sprintf("El error relativo de integración por el método de Monte Carlo es de %0.7g",eRel))