using Plots, LaTeXStrings

function mDot(xi::Number,mi::Number,K::Number,c1::Number,c2::Number,p::Number,q::Number,delta::Number,r::Number)
    y=(delta-r)*mi + 2*r*mi*xi/K - (p-mi)^2*q^2*xi/c2 + (p-mi)*q*c1/c2
    return y
end

function xDot(xi::Number,mi::Number,K::Number,c1::Number,c2::Number,p::Number,q::Number,delta::Number,r::Number)
    y=r*xi*(1-xi/K) - (p-mi)*(q*xi)^2/c2 + xi*q*c1/c2
    return y
end

function fishery(m0::Number,K::Number,c1::Number,c2::Number,p::Number,q::Number,delta::Number,r::Number,n::Int,x::AbstractVector,m::AbstractVector)
    h=1/n
    m[1]=m0
    for i=1:n
        x[i+1]=x[i]+ h*xDot(x[i],m[i],K,c1,c2,p,q,delta,r)
        m[i+1]=m[i]+ h*mDot(x[i],m[i],K,c1,c2,p,q,delta,r)
    end
end

function fisheryStochasticApprox(K::Number,c1::Number,c2::Number,p::Number,q::Number,delta::Number,r::Number,n::Int,tol::AbstractFloat)
    x=zeros(Float64,n+1)
    x[1]=K
    m=copy(x)
    m0=0.0
    i=1
    finalCondition=0.0
    notSolved=true
    while notSolved && i<1000
        fishery(m0,K,c1,c2,p,q,delta,r,n,x,m)
        finalCondition=m[n+1]
        notSolved=abs(finalCondition)>tol
        m0-=(finalCondition)/i
        i+=1
    end
    return x,m
end

function fisheryShooting(K::Number,c1::Number,c2::Number,p::Number,q::Number,delta::Number,r::Number,n::Int,tol::AbstractFloat)
    x=zeros(Float64,n+1)
    x[1]=K
    x2=copy(x)
    m=copy(x)
    m2=copy(x)
    m0=0.3
    i=1
    dm=1.0
    finalCondition=0.0
    notSolved=true
    while notSolved && i<1000
        fishery(m0,K,c1,c2,p,q,delta,r,n,x,m)
        fishery(m0+1/n,K,c1,c2,p,q,delta,r,n,x2,m2)
        finalCondition=m[n+1]
        dm=(m2[n+1]-m[n+1])*n
        notSolved=abs(finalCondition)>tol
        m0=m0-finalCondition/dm
        i+=1
    end
    return x,m
end

r=0.71
delta=0.12
q=0.0001
K=10^6
p=0.5


c1=0.001
c2=0.001
n=1000
tol=1e-10
x,m=fisheryStochasticApprox(K,c1,c2,p,q,delta,r,n,tol)

E=((p.-m).*q.*x .- c1)/c2
t=LinRange(0,1,n+1)

plt1=plot(t,E,color_palette=:lake,markersize=2,label=L"E(x,t)",lw = 1.5,size=(800,600),legend=:bottomright,ylims=(7100,18000),xlims=(0,1));
savefig(plt1,"Ehamiltoniano.png") #Guardamos la gráfica