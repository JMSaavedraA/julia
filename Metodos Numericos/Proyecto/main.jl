include("Solvers.jl")
include("metodosTridiagonales.jl")
include("sitnikov.jl")
using Plots
using LaTeXStrings
using Printf

#Parámetros
p=1
c=2
exc=0.20
n=1000000
T=10
t=LinRange(0,T,n+1)
z0=1.25
z0dot=0.0
tol=1e-8
maxIter=1000
h=t[2]-t[1]
theta=zeros(Float64,n+1)
#Resolvemos θ'=4c(1+ϵcos(θ))^2/p^2
@time eulerTheta(n,theta,h,c,exc,p)
#Calculamos r=p/(2(1+ϵcos(θ)))
R=(p/2)./(1 .+ exc*cos.(theta))
#Resolvemos z"=-z/(z^2+r^2)^()
@time Z=sitnikov(z0,z0dot,R,n,tol,maxIter,h)

plt1=plot(t,Z,label = L"z",legend=:topleft,size=(800,600)); #Graficamos
savefig(plt1,@sprintf("sitnikov%0.2f.png",z0)) #Guardamos la gráfica

plt2=plot(t,theta,label = L"\theta",legend=:topleft,size=(800,600)); #Graficamos
savefig(plt2,@sprintf("theta%0.2f.png",z0)) #Guardamos la gráfica

plt3=plot(t,R,label = L"r",legend=:topleft,size=(800,600)); #Graficamos
savefig(plt3,@sprintf("r%0.2f.png",z0)) #Guardamos la gráfica