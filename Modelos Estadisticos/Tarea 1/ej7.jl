using Plots, Distributions, SpecialFunctions, LaTeXStrings

t=LinRange(0.05,10,1000)
actual=2*(1 .- cdf.(Normal(),t))
function MarkovBound(t::Number, k::Int)
    2^(k/2)*gamma((k+1)/2)/(sqrt(pi)*t^k)
end
markov1=MarkovBound.(t,1)
markov2=MarkovBound.(t,2)
markov3=MarkovBound.(t,3)
markov4=MarkovBound.(t,4)
markov5=MarkovBound.(t,5)

mills=sqrt(2/pi)*exp.(-t.^2 ./ 2) ./ t

plt1=plot(t,actual,ylim=(0, 1),size=(800,600),color_palette=:seaborn_bright6,label=L"\mathbb{P}(|Z|>t)",lw=1.5);
plot!(t,mills,color_palette=:seaborn_bright6, label="Mill's bound",lw=1.5);
plot!(t,[markov1,markov2,markov3,markov4,markov5],color_palette=:BuPu_7, label=[L"k=1" L"k=2" L"k=3" L"k=4" L"k=5"],lw=1.5);

savefig(plt1,"ej7.png")