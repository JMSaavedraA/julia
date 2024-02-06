using Printf
include("gaussLegendre.jl")
int(x) = floor(Int, x)

#Definimos la función f_1(x)
function f1(x::Number)
    y=(1+x^2)^(-2)
    return y
end

#Límites de integración
a=-1.0
b=1.0
m=100 #h=2/m cambiar a m=10, m=100, m=1000 para reproducir las tablas del pdf.
#Inicializamos
I=zeros(Float64,3)
eAbs=copy(I)
eRel=copy(I)

println(@sprintf("Se considera la cuadratura de Gauss-Legendre con m=%i intervalos",m))
for i=1:3
    #Calculamos la integral de f_1(x) en D
    I[i]=cuadraturaGaussiana(f1,a,b,i,m)
    #Calculamos error absoluto y relativo
    eAbs[i]=abs(I[i]-(pi+2)/4)
    eRel[i]=eAbs[i]*4/(pi+2)
    println(@sprintf("Para n=%i, es decir, %i puntos por intervalo",i-1,2*i-1))
    println(@sprintf("El valor aproximado de la integral es de %0.7g",I[i]))
    println(@sprintf("El error absoluto de integración es de %0.7g",eAbs[i]))
    println(@sprintf("El error relativo de integración es de %0.7g",eRel[i]))
end