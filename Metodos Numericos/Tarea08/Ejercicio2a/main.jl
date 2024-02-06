include("Solvers.jl")
include("interpolacion.jl")
using Random

Random.seed!(2022); #Tomamos una semilla aleatoria para poder replicar los resultados

m=1000
z=sort(2*rand(m).-1) #z serán los 1000 puntos en que se evalua el polinomio resultante y la función original
f1=(25*z.^2 .+ 1).^-1
f2=abs.(z)-z/2 - z.^2

n=3 #Número de puntos donde se interpola el polinomio
x=sort(2*rand(n).-1) #Los puntos sobre los que se interpola el polinomio
y1=(25*x.^2 .+ 1).^-1 
y2=abs.(x)-x/2 - x.^2


p1=interpolarPolinomios(z,x,y1,m,n) #Interpolamos en n puntos y se evalua en z
p2=interpolarPolinomios(z,x,y2,m,n)

ecm1=norma2(p1-f1,m)/m; #Cálculo del Error Cuadrático Medio
ecm2=norma2(p2-f2,m)/m;
println(@sprintf("El error cuadrático medio de la aproximación de f1 interpolando por %i puntos es de %g",n,ecm1))
println(@sprintf("El error cuadrático medio de la aproximación de f2 interpolando por %i puntos es de %g",n,ecm2))
guardarVector(p1,m,@sprintf("p1-%i.txt",n))
guardarVector(p2,m,@sprintf("p2-%i.txt",n))

n=5 #Número de puntos donde se interpola el polinomio
x=sort(2*rand(n).-1) #Los puntos sobre los que se interpola el polinomio
y1=(25*x.^2 .+ 1).^-1
y2=abs.(x)-x/2 - x.^2

p1=interpolarPolinomios(z,x,y1,m,n)
p2=interpolarPolinomios(z,x,y2,m,n)

ecm1=norma2(p1-f1,m)/m; #Cálculo del Error Cuadrático Medio
ecm2=norma2(p2-f2,m)/m;
println(@sprintf("El error cuadrático medio de la aproximación de f1 interpolando por %i puntos es de %g",n,ecm1))
println(@sprintf("El error cuadrático medio de la aproximación de f2 interpolando por %i puntos es de %g",n,ecm2))
guardarVector(p1,m,@sprintf("p1-%i.txt",n))
guardarVector(p2,m,@sprintf("p2-%i.txt",n))

n=8 #Número de puntos donde se interpola el polinomio
x=sort(2*rand(n).-1) #Los puntos sobre los que se interpola el polinomio
y1=(25*x.^2 .+ 1).^-1
y2=abs.(x)-x/2 - x.^2

p1=interpolarPolinomios(z,x,y1,m,n)
p2=interpolarPolinomios(z,x,y2,m,n)

ecm1=norma2(p1-f1,m)/m; #Cálculo del Error Cuadrático Medio
ecm2=norma2(p2-f2,m)/m;
println(@sprintf("El error cuadrático medio de la aproximación de f1 interpolando por %i puntos es de %g",n,ecm1))
println(@sprintf("El error cuadrático medio de la aproximación de f2 interpolando por %i puntos es de %g",n,ecm2))
guardarVector(p1,m,@sprintf("p1-%i.txt",n))
guardarVector(p2,m,@sprintf("p2-%i.txt",n))