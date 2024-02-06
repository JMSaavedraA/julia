include("Splines.jl")
include("interpolacion.jl")
using Printf
using Plots
using LaTeXStrings
using Random

Random.seed!(2022);

m=1000 #Consideramos mil puntos uniformes en (0,50) para medir el error
z=LinRange(0,50,m)
f1=z.+z.*sin.(z/2)/3
yMax=maximum(f1)
yMin=minimum(f1)
ddI=4/3
ddD=1+(cos(25)-12.5*sin(25))/3


n=30 #Tomamos 30 puntos aleatorios uniformes (0,50)
x=sort(50*rand(n))
x[1]=0 #Tomamos los extremos del intervalo
x[n]=50
y1=x.+x.*sin.(x/2)/3

pI1=interpolarPolinomios(z,x,y1,m,n) #Calculamos la interpolación por polinomio de grado n-1 por los puntos (x,y) evaluada en z
ecmI1=norma2(pI1-f1,m)/m; #Cálculo del Error Cuadrático Medio para el polinomio de grado n-1
pL1=polinomioLagrange(z,x,y1,m,n) #Calculamos la interpolación por polinomio de Lagrange por los puntos (x,y) evaluada en z
ecmL1=norma2(pL1-f1,m)/m; #Cálculo del Error Cuadrático Medio para el polinomio de Lagrange
pN1=polinomioNewton(z,x,y1,m,n) #Calculamos la interpolación por polinomio de Newton por los puntos (x,y) evaluada en z
ecmN1=norma2(pN1-f1,m)/m; #Cálculo del Error Cuadrático Medio para el polinomio de Newton
pc01=splineLineal(z,x,y1,m,n) #Calculamos la interpolación por Splines de continuidad 0 por los puntos (x,y) evaluada en z
ecmc01=norma2(pc01-f1,m)/m; #Cálculo del Error Cuadrático Medio de los Splines de continuidad 0
pc11=splineCuadraticoIzquierda(z,x,y1,m,n,1.0) #Calculamos la interpolación por Splines de continuidad 1 por los puntos (x,y) evaluada en z
ecmc11=norma2(pc11-f1,m)/m; #Cálculo del Error Cuadrático Medio de los Splines de continuidad 1
pc21=splineCubico(z,x,y1,m,n,ddI,ddD) #Calculamos la interpolación por Splines de continuidad 2 por los puntos (x,y) evaluada en z
ecmc21=norma2(pc21-f1,m)/m; #Cálculo del Error Cuadrático Medio de los Splines de continuidad 2

plt1=plot(z,[f1,pI1,pL1,pN1,pc01,pc11,pc21],label = [L"x + \frac{x\sin(\frac{x}{2})}{3}" "polinomio interpolado" "polinomio de Lagrange" "polinomio de Newton" "spline c=0" "spline c=1" "spline c=2"],legend=:topleft,ylims=(yMin,yMax),size=(800,600));scatter!(x,y1,label="puntos interpolados")
savefig(plt1,@sprintf("comparacion-%i.png",n)) #Graficamos

println(@sprintf("El error cuadrático medio de la aproximación de f1 interpolando con polinomio de grado %i es de %g",n-1,ecmI1))
guardarVector(pI1,m,@sprintf("pInterpolado-%i.txt",n))
println(@sprintf("El error cuadrático medio de la aproximación de f1 interpolando con polinomio de Lagrange por %i puntos es de %g",n,ecmL1))
guardarVector(pL1,m,@sprintf("pLagrange-%i.txt",n))
println(@sprintf("El error cuadrático medio de la aproximación de f1 interpolando con polinomio de Newton por %i puntos es de %g",n,ecmN1))
guardarVector(pN1,m,@sprintf("pNewton-%i.txt",n))
println(@sprintf("El error cuadrático medio de la aproximación de f1 interpolando con splines lineales por %i puntos es de %g",n,ecmc01))
guardarVector(pc01,m,@sprintf("pLineal-%i.txt",n))
println(@sprintf("El error cuadrático medio de la aproximación de f1 interpolando con splines cuadráticos por %i puntos es de %g",n,ecmc11))
guardarVector(pc11,m,@sprintf("pCuadratico-%i.txt",n))
println(@sprintf("El error cuadrático medio de la aproximación de f1 interpolando con splines cúbicos por %i puntos es de %g",n,ecmc21))
guardarVector(pc21,m,@sprintf("pCubicos-%i.txt",n))