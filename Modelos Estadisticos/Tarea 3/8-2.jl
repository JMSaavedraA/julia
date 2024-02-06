using Random, Statistics, Distributions, StatsBase

Random.seed!(2023)
m=1000;
n=50;
B=1000;
b=zeros(Float64,B,m)
lB=zeros(Float64,3)
uB=zeros(Float64,3)
te = trues(3,m)
za=quantile(Normal(),0.975)

actualSkew = (exp(1)+2)*sqrt(exp(1)-1)

for i=1:m
    x=rand(LogNormal(),n);
    for j=1:B
        y=sample(x,n);
        b[j,i]=skewness(y)
    end
    se=std(b[:,i])
    est=skewness(x)
    lB[1]=est - za*se
    uB[1]=est + za*se
    lB[2]=2*est - quantile(b[:,i],0.975)
    uB[2]=2*est + quantile(b[:,i],0.025)
    lB[3]=quantile(b[:,i],0.025)
    uB[3]=quantile(b[:,i],0.975)
    te[:,i] = lB.<=actualSkew.<=uB
end

correctGuess = sum(te,dims=2)./m