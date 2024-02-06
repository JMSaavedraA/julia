using Random
using Printf
include("Integracion.jl")

Random.seed!(2022)

#Definimos la función f_1(x)
function f1(x::AbstractFloat)
    y=(1+x^2)^(-2)
    return y
end

#Límites de integración
a=-1.0
b=1.0
#valores superiores e inferiores de f_1(x)
fMax=1.0
fMin=0.25
#Calculamos la integral de f_1(x) en D
I=integralMonteCarlo(f1,a,b,10000,1000,1e-15,fMax,fMin)
#Error absoluto y relativo
eAbs=abs(I-(pi+2)/4)
eRel=eAbs*4/(pi+2)

println(@sprintf("El valor aproximado de la integral es de %0.7g",I))
println(@sprintf("El error absoluto de integración por el método de Monte Carlo es de %0.7g",eAbs))
println(@sprintf("El error relativo de integración por el método de Monte Carlo es de %0.7g",eRel))