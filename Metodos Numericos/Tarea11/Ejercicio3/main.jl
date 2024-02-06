using Printf
include("extrapolacion.jl")
int(x) = floor(Int, x)

#Definimos la función f_1(x)
function f1(x::AbstractFloat)
    y=(1+x^2)^(-2)
    return y
end

#Límites de integración
a=-1.0
b=1.0
#Calculamos la integral de f_1(x) en D
m=10 #m=5, 8 y 10 en el PDF
I=integralRomberg(f1,a,b,m)

#Error absoluto y relativo
eAbs=abs(I-(pi+2)/4)
eRel=eAbs*4/(pi+2)

println(@sprintf("Se considera la extrapolación de Richardson con m=%i",m))
println(@sprintf("El valor aproximado de la integral es de %0.7g", I))
println(@sprintf("El error absoluto de integración es de %0.7g", eAbs))
println(@sprintf("El error relativo de integración es de %0.7g", eRel))