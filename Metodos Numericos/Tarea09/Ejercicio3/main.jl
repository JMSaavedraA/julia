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

n=12 #Primero tomamos 10 puntos aleatorios uniformes (0,50)
x=sort(50*rand(n))
x[1]=0 #Tomamos los extremos del intervalo
x[n]=50
y1=x.+x.*sin.(x/2)/3

p1=splineCubico(z,x,y1,m,n,ddI,ddD) #Calculamos la interpolación por Splines por los puntos (x,y) evaluada en z

plt1=plot(z,[f1,p1],label = [L"x + \frac{x\sin(\frac{x}{2})}{3}" "spline c=2"],legend=:topleft,ylims=(yMin,yMax),size=(800,600)); #Graficamos
savefig(plt1,@sprintf("pSplineCubico1-%i.png",n)) #Guardamos la gráfica
ecm1=norma2(p1-f1,m)/m; #Cálculo del Error Cuadrático Medio
println(@sprintf("El error cuadrático medio de la aproximación de f1 interpolando con splines cúbicos por %i puntos es de %g",n,ecm1))
guardarVector(p1,m,@sprintf("p1-%i.txt",n)) #Guardamos el vector


n=102 #Ahora aumentamos a 100 puntos aleatorios uniformes (0,50)
x=sort(50*rand(n))
x[1]=0 #Tomamos los extremos del intervalo
x[n]=50
y1=x.+x.*sin.(x/2)/3

p1=splineCubico(z,x,y1,m,n,ddI,ddD) #Calculamos la interpolación por Splines por los puntos (x,y) evaluada en z

plt1=plot(z,[f1,p1],label = [L"x + \frac{x\sin(\frac{x}{2})}{3}" "spline c=2"],legend=:topleft,ylims=(yMin,yMax),size=(800,600)); #Graficamos
savefig(plt1,@sprintf("pSplineCubico1-%i.png",n)) #Guardamos la gráfica
ecm1=norma2(p1-f1,m)/m; #Cálculo del Error Cuadrático Medio
println(@sprintf("El error cuadrático medio de la aproximación de f1 interpolando con splines cúbicos por %i puntos es de %g",n,ecm1))
guardarVector(p1,m,@sprintf("p1-%i.txt",n)) #Guardamos el vector


n=1002 #Finalmente tomamos 1000 puntos aleatorios uniformes
x=sort(50*rand(n))
x[1]=0 #Tomamos los extremos del intervalo
x[n]=50
y1=x.+x.*sin.(x/2)/3

p1=splineCubico(z,x,y1,m,n,ddI,ddD) #Calculamos la interpolación por Splines por los puntos (x,y) evaluada en z

plt1=plot(z,[f1,p1],label = [L"x + \frac{x\sin(\frac{x}{2})}{3}" "spline c=2"],legend=:topleft,ylims=(yMin,yMax),size=(800,600)); #Graficamos
savefig(plt1,@sprintf("pSplineCubico1-%i.png",n)) #Guardamos la gráfica
ecm1=norma2(p1-f1,m)/m; #Cálculo del Error Cuadrático Medio
println(@sprintf("El error cuadrático medio de la aproximación de f1 interpolando con splines cúbicos por %i puntos es de %g",n,ecm1))
guardarVector(p1,m,@sprintf("p1-%i.txt",n)) #Guardamos el vector