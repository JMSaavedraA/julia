using Plots, LaTeXStrings

function EDot(xi::Number,Ei::Number,K::Number,c1::Number,c2::Number,p::Number,q::Number,delta::Number,r::Number)
    y=r*Ei*xi/K + delta*Ei + p*q*(r-delta)*xi/c2 - 2*r*p*q*xi^2/(c2*K) + c1*(delta+r*xi/K)/c2
    return y
end

function xDot(xi::Number,Ei::Number,K::Number,c1::Number,c2::Number,p::Number,q::Number,delta::Number,r::Number)
    y=r*xi*(1-xi/K) - q*Ei*xi
    return y
end

function fishery(E0::Number,K::Number,c1::Number,c2::Number,p::Number,q::Number,delta::Number,r::Number,n::Int,x::AbstractVector,E::AbstractVector)
    h=1/n
    E[1]=E0
    for i=1:n
        x[i+1]=x[i]+ h*xDot(x[i],E[i],K,c1,c2,p,q,delta,r)
        E[i+1]=E[i]+ h*EDot(x[i],E[i],K,c1,c2,p,q,delta,r)
    end
end

function fisheryShooting(K::Number,c1::Number,c2::Number,p::Number,q::Number,delta::Number,r::Number,n::Int,tol::AbstractFloat)
    x=zeros(Float64,n+1)
    x[1]=K
    x2=copy(x)
    E=copy(x)
    E2=copy(x)
    E0=17500
    i=1
    dE=1.0
    finalCondition=0.0
    notSolved=true
    while notSolved && i<1000
        fishery(E0,K,c1,c2,p,q,delta,r,n,x,E)
        fishery(E0+1/n,K,c1,c2,p,q,delta,r,n,x2,E2)
        finalCondition=E[n+1]-(p*q*x[n+1]-c1)/c2
        #dE=(E2[n+1]-E[n+1])*n
        dE=(E2[n+1]-E[n+1] + p*q*(x[n+1]-x2[n+1]))*n
        notSolved=abs(finalCondition)>tol
        E0=E0-finalCondition/dE
        i+=1
    end
    return x,E
end

r=0.71
delta=0.12
q=0.0001
K=10^6
p=0.5

c1=0.001
c2=0.001
n=1000
tol=1e-9
x,E=fisheryShooting(K,c1,c2,p,q,delta,r,n,tol)
t=LinRange(0,1,n+1)

plt1=plot(t,E,color_palette=:lake,markersize=2,label=L"E(x,t)",lw = 1.5,size=(800,600),legend=:bottomright,ylims=(7100,18000),xlims=(0,1));
savefig(plt1,"Evariacional.png") #Guardamos la gráfica