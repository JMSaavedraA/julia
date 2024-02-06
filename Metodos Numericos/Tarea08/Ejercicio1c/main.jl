include("Solvers.jl")
include("interpolacion.jl")
using Random

Random.seed!(2022); #Tomamos una semilla aleatoria para poder replicar los resultados

function g1(x::AbstractFloat) #Proponemos 3 funciones para aproximar f1
    g=[1,exp(-x^2),cospi(x)]
    return g
end

function g2(x::AbstractFloat) #Proponemos 3 funciones para aproximar f2
    g=[x,exp(-x^2),cospi(2*x)]
    return g
end

n=3
m=5
z=sort(2*rand(m).-1)
y1=(25*z.^2 .+ 1).^-1
y2=abs.(z)-z/2 - z.^2


p1=minimosCuadrados(g1,z,y1,m,n) #Hacemos minimos cuadrados para aproximación con g1 a f1 con 5 puntos
p2=minimosCuadrados(g2,z,y2,m,n) #Hacemos minimos cuadrados para aproximación con g2 a f2 con 5 puntos

ecm1=norma2(p1-y1,m)/m; #Cálculo del Error Cuadrático Medio
ecm2=norma2(p2-y2,m)/m;
println(@sprintf("El error cuadrático medio de la aproximación de f1 para %i puntos es de %g",m,ecm1))
println(@sprintf("El error cuadrático medio de la aproximación de f2 para %i puntos es de %g",m,ecm2))
guardarVector(p1,m,@sprintf("p1-%i.txt",m))
guardarVector(p2,m,@sprintf("p2-%i.txt",m))

m=100
z=sort(2*rand(m).-1)
y1=(25*z.^2 .+ 1).^-1
y2=abs.(z)-z/2 - z.^2

p1=minimosCuadrados(g1,z,y1,m,n) #Hacemos minimos cuadrados para aproximación con g1 a f1 con 100 puntos
p2=minimosCuadrados(g2,z,y2,m,n) #Hacemos minimos cuadrados para aproximación con g2 a f2 con 100 puntos

ecm1=norma2(p1-y1,m)/m; #Cálculo del Error Cuadrático Medio
ecm2=norma2(p2-y2,m)/m;
println(@sprintf("El error cuadrático medio de la aproximación de f1 para %i puntos es de %g",m,ecm1))
println(@sprintf("El error cuadrático medio de la aproximación de f2 para %i puntos es de %g",m,ecm2))
guardarVector(p1,m,@sprintf("p1-%i.txt",m))
guardarVector(p2,m,@sprintf("p2-%i.txt",m))

m=1000
z=sort(2*rand(m).-1)
y1=(25*z.^2 .+ 1).^-1
y2=abs.(z)-z/2 - z.^2

p1=minimosCuadrados(g1,z,y1,m,n) #Hacemos minimos cuadrados para aproximación con g1 a f1 con 1000 puntos
p2=minimosCuadrados(g2,z,y2,m,n) #Hacemos minimos cuadrados para aproximación con g1 a f1 con 1000 puntos

ecm1=norma2(p1-y1,m)/m; #Cálculo del Error Cuadrático Medio
ecm2=norma2(p2-y2,m)/m;
println(@sprintf("El error cuadrático medio de la aproximación de f1 para %i puntos es de %g",m,ecm1))
println(@sprintf("El error cuadrático medio de la aproximación de f2 para %i puntos es de %g",m,ecm2))
guardarVector(p1,m,@sprintf("p1-%i.txt",m))
guardarVector(p2,m,@sprintf("p2-%i.txt",m))