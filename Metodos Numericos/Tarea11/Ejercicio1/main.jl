using Printf
include("newtonCotes.jl")
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
m=10 #h=2/m cambiar a m=12, m=120 para reproducir las tablas del pdf.
I=zeros(Float64,4)
I[1]=integralTrapecio(f1,a,b,m)
I[2]=integralSimpson(f1,a,b,m)
I[3]=integralTresOctavos(f1,a,b,m)
I[4]=integralMilne(f1,a,b,m)
#Error absoluto y relativo
eAbs=abs.(I.-(pi+2)/4)
eRel=eAbs*4/(pi+2)

println(@sprintf("Se considera el método de Newton-Cotes con m=%i",m))
for i=1:4
    println(@sprintf("Para Newton-Cotes con n=%i",i))
    println(@sprintf("El valor aproximado de la integral es de %0.7g",I[i]))
    println(@sprintf("El error absoluto de integración es de %0.7g",eAbs[i]))
    println(@sprintf("El error relativo de integración es de %0.7g",eRel[i]))
end