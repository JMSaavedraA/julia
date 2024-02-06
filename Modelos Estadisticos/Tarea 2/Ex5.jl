using PlotlyJS, Distributions, Statistics, Random, StatsBase

Random.seed!(2023)
n=100
m=1000
k=10000
t=LinRange(-5,5,k)
y=cdf.(Cauchy(),t)
x=rand(Cauchy(),n,m)
epsilon=sqrt(log(2/0.95)/(2n))
percent=zeros(m)
for i=1:m
    estF=ecdf(x[:,i])
    L(x)=maximum([estF(x)-epsilon,0])
    U(x)=minimum([estF(x)+epsilon,1])
    percent[i]=sum(L.(t) .< y .< U.(t))/k
end

plot(box(y=percent, quartilemethod="linear", name="Correct Prediction Ratio"),config=PlotConfig(responsive=false))