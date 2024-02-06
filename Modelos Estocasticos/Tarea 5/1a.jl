using Random
using Plots
using LaTeXStrings

Random.seed!(2022)

c=100
lambda=5
Nc=121
a=randexp(Nc)
T=cumsum(a)
t=LinRange(0,c,10000)
j=1
d=t[1]
N=zeros(Int64,10000).+Nc
for i=1:Nc
    r=T[i]
    while d<r && j<10000
        N[j]=i-1
        global j+=1
        global d=t[j]
    end
end

X=N-t*Nc/c

plt1=scatter(T,1:121,color_palette=:lake,markersize=2,label=L"N(T_i)");plot!(t,N,lw = 1.5, label=L"N(t)",size=(1600,1200),legend=:bottomright);
savefig(plt1,"N.png") #Guardamos la gráfica

plt2=plot(t,X,lw = 1,color_palette=:dracula, label=L"X(t)",size=(1600,1200),legend=:bottomright);
savefig(plt2,"X.png") #Guardamos la gráfica