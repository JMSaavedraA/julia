using DelimitedFiles, Statistics, Distributions, PlotlyJS
A=readdlm("oldFaithful.txt",'\t')
waitTimes=A[:,2]

xBar=mean(waitTimes)
se=std(waitTimes,mean=xBar)
zAlpha=quantile(Normal(),0.95)

lBound=xBar - zAlpha*se
uBound=xBar + zAlpha*se

xMedian=median(waitTimes)

sum(lBound.<waitTimes .< uBound)/length(waitTimes)*100