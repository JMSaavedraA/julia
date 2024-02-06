include("elementoFinito.jl")
include("Solvers.jl")
using Printf
using Plots
using LaTeXStrings
using Random

Random.seed!(2022);
m=2001
k=501
n=101
fuzzy=1.0
x=LinRange(-5,5,k)
y=sinpi.(x) + fuzzy*rand(k) .- fuzzy/2 #y_k es sin(pi*x_k) + u_k, donde u_k es aleatorio
z=LinRange(-5,5,m)
fz=sinpi.(z)
w=LinRange(-5,5,n)

#Primero, consideramos λ=0.5

lambda=0.5
fz1=elementoFinitoLineal(z,x,y,w,lambda)
ecm1=norma2(fz-fz1,m)/m; #Cálculo del Error Cuadrático Medio
println(@sprintf("El error cuadrático medio de la aproximación por elementos finitos con λ=%1.1f es de %g",lambda,ecm1))
guardarVector(fz1,m,@sprintf("elementoFinito-%1.1f.txt",lambda))

#Ahora, consideramos λ=1.5

lambda=1.5
fz2=elementoFinitoLineal(z,x,y,w,lambda)
ecm2=norma2(fz-fz2,m)/m; #Cálculo del Error Cuadrático Medio
println(@sprintf("El error cuadrático medio de la aproximación por elementos finitos con λ=%1.1f es de %g",lambda,ecm2))
guardarVector(fz2,m,@sprintf("elementoFinito-%1.1f.txt",lambda))

#Finalmente, consideramos λ=3.5

lambda=3.5
fz3=elementoFinitoLineal(z,x,y,w,lambda)
ecm3=norma2(fz-fz3,m)/m; #Cálculo del Error Cuadrático Medio
println(@sprintf("El error cuadrático medio de la aproximación por elementos finitos con λ=%1.1f es de %g",lambda,ecm3))
guardarVector(fz3,m,@sprintf("elementoFinito-%1.1f.txt",lambda))

plt1=scatter(x,y,seriescolor=:Blues,label=L"(x_k,y_k)"); plot!(z,[fz fz1 fz2 fz3],lw = 3,color_palette=:seaborn_bright6, label=[L"\sin(\pi x)" L"\lambda=0.5" L"\lambda=1.5" L"\lambda=3.5"],size=(800,600),legend=:bottomright);
savefig(plt1,"elementoFinitoComparacion.png") #Guardamos la gráfica