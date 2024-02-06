include("Solvers.jl")
include("interpolacion.jl")
using Random

Random.seed!(2022); #Tomamos una semilla aleatoria para poder replicar los resultados

n=1 #Hacemos minimos cuadrados para aproximación lineal de f1 y f2 con 5 puntos
m=5
z=sort(2*rand(m).-1)
y1=(25*z.^2 .+ 1).^-1
y2=abs.(z)-z/2 - z.^2

p1=minimosCuadradosPolinomios(z,y1,m,n)
p2=minimosCuadradosPolinomios(z,y2,m,n)

ecm1=norma2(p1-y1,m)/m; #Cálculo del Error Cuadrático Medio
ecm2=norma2(p2-y2,m)/m;
println(@sprintf("El error cuadrático medio de la aproximación de f1 para %i puntos es de %g",m,ecm1))
println(@sprintf("El error cuadrático medio de la aproximación de f2 para %i puntos es de %g",m,ecm2))
guardarVector(p1,m,@sprintf("p1-%i.txt",m))
guardarVector(p2,m,@sprintf("p2-%i.txt",m))

m=100 #Hacemos minimos cuadrados para aproximación lineal de f1 y f2 con 100 puntos
z=sort(2*rand(m).-1)
y1=(25*z.^2 .+ 1).^-1
y2=abs.(z)-z/2 - z.^2

p1=minimosCuadradosPolinomios(z,y1,m,n)
p2=minimosCuadradosPolinomios(z,y2,m,n)

ecm1=norma2(p1-y1,m)/m; #Cálculo del Error Cuadrático Medio
ecm2=norma2(p2-y2,m)/m;
println(@sprintf("El error cuadrático medio de la aproximación de f1 para %i puntos es de %g",m,ecm1))
println(@sprintf("El error cuadrático medio de la aproximación de f2 para %i puntos es de %g",m,ecm2))
guardarVector(p1,m,@sprintf("p1-%i.txt",m))
guardarVector(p2,m,@sprintf("p2-%i.txt",m))

m=1000 #Hacemos minimos cuadrados para aproximación lineal de f1 y f2 con 1000 puntos
z=sort(2*rand(m).-1)
y1=(25*z.^2 .+ 1).^-1
y2=abs.(z)-z/2 - z.^2

p1=minimosCuadradosPolinomios(z,y1,m,n)
p2=minimosCuadradosPolinomios(z,y2,m,n)

ecm1=norma2(p1-y1,m)/m; #Cálculo del Error Cuadrático Medio
ecm2=norma2(p2-y2,m)/m;
println(@sprintf("El error cuadrático medio de la aproximación de f1 para %i puntos es de %g",m,ecm1))
println(@sprintf("El error cuadrático medio de la aproximación de f2 para %i puntos es de %g",m,ecm2))
guardarVector(p1,m,@sprintf("p1-%i.txt",m))
guardarVector(p2,m,@sprintf("p2-%i.txt",m))