using Distributions
using Random
using Plots
using LaTeXStrings
using Printf
using KernelDensity

Random.seed!(2022)

k=100
d=Poisson(4)
x=rand(d,k)
m=mean(x)
function logVerosimilitud(l::Number)
    y=-l*k + (k*m)*log(l) - sum(log.(factorial.(x)))
    return y
end
function noisyLikelihood(l::Number)
    y=logVerosimilitud(l)+rand(Normal(0,0.5))
    return y
end
function errorTruncamiento(l::Number)
    y=logVerosimilitud(l) - noisyLikelihood(l)
    return y
end

n=3000
N=10000
tol=1e-5
r=100
l=zeros(Float64,n+1)
estimador=zeros(Float64,r)
l[1]=2.0
h=1e-3
for j=1:r
    for i=1:n
        dl=(noisyLikelihood(l[i]+h)-noisyLikelihood(l[i]-h))/(2*h)
        if abs(dl)<tol
            break
        end
        l[i+1]=minimum([maximum([l[i]+dl/(i*log(N)),1]),10])
    end
    estimador[j]=l[n+1]
end

media=mean(estimador)
desv=sqrt(var(estimador))
estimador=(estimador.-media)./desv
U=kde(estimador)

plt1=histogram(estimador,normalize = :probability,bins=:sturges,label=L"\hat{x}",size=(1600,1200),legend=:topright,palette=:Reds_3); plot!(U.x,pdf.(Normal(mean(estimador), sqrt(var(estimador))),U.x),lw=5,label="Densidad Normal");plot!(U.x,U.density,lw=5,label="Estimación por Kernel Normal");
savefig(plt1,@sprintf("histograma3-%i.png",n)) #Guardamos la gráfica