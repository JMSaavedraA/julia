using Random, Statistics, Distributions, StatsBase, Plots, LaTeXStrings

Random.seed!(2023)
n=100;
B=1000;
b=zeros(Float64,B)
lB=zeros(Float64,3)
uB=zeros(Float64,3)
za=quantile(Normal(),0.975)


x=rand(Normal(5,1),n);
for j=1:B
    y=sample(x,n);
    b[j]=exp(mean(y))
end
se=std(b)
est=exp(mean(x))
lB[1]=est - za*se
uB[1]=est + za*se
lB[2]=2*est - quantile(b,0.975)
uB[2]=2*est + quantile(b,0.025)
lB[3]=quantile(b,0.025)
uB[3]=quantile(b,0.975)
t=LinRange(0,maximum(uB),B)

plt1=histogram(b,normalize = :pdf,label=L"\hat{\theta}",size=(1600,1200),legend=:topright,palette=:Oranges_3);
plot!(t,pdf.(LogNormal(5, 1/sqrt(n)),t),lw=5,label="Densidad LogNormal")
savefig(plt1,"8-6.png")