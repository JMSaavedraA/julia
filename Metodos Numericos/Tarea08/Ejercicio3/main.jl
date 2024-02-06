include("interpolacion.jl")
using Plots
using LaTeXStrings
using Random

Random.seed!(2022);

m=1000
z=sort(2*rand(m).-1)
f1=(25*z.^2 .+ 1).^-1
f2=abs.(z)-z/2 - z.^2

n=3
x=sort(2*rand(n).-1)
y1=(25*x.^2 .+ 1).^-1
y2=abs.(x)-x/2 - x.^2


pI1=interpolarPolinomios(z,x,y1,m,n)
pI2=interpolarPolinomios(z,x,y2,m,n)
pL1=polinomioLagrange(z,x,y1,m,n)
pL2=polinomioLagrange(z,x,y2,m,n)
pN1=polinomioNewton(z,x,y1,m,n)
pN2=polinomioNewton(z,x,y2,m,n)

plt1=plot(z,[f1,pI1,pL1,pN1],label = [L"\frac{1}{1+25x^2}" "polinomio interpolado" "polinomio de Lagrange" "polinomio de Newton"],legend=:topleft); scatter!(x,y1,label="puntos Interpolados");
plt2=plot(z,[f2,pI2,pL2,pN2],label = [L"|x| - \frac{x}{2} - x^2" "polinomio interpolado" "polinomio de Lagrange" "polinomio de Newton"],legend=:bottomleft); scatter!(x,y2,label="puntos Interpolados");
savefig(plt1,@sprintf("pInterpolacion1-%i.png",n))
savefig(plt2,@sprintf("pInterpolacion2-%i.png",n))

n=5
x=sort(2*rand(n).-1)
y1=(25*x.^2 .+ 1).^-1
y2=abs.(x)-x/2 - x.^2

pI1=interpolarPolinomios(z,x,y1,m,n)
pI2=interpolarPolinomios(z,x,y2,m,n)
pL1=polinomioLagrange(z,x,y1,m,n)
pL2=polinomioLagrange(z,x,y2,m,n)
pN1=polinomioNewton(z,x,y1,m,n)
pN2=polinomioNewton(z,x,y2,m,n)

plt1=plot(z,[f1,pI1,pL1,pN1],label = [L"\frac{1}{1+25x^2}" "polinomio interpolado" "polinomio de Lagrange" "polinomio de Newton"],legend=:topleft); scatter!(x,y1,label="puntos Interpolados");
plt2=plot(z,[f2,pI2,pL2,pN2],label = [L"|x| - \frac{x}{2} - x^2" "polinomio interpolado" "polinomio de Lagrange" "polinomio de Newton"],legend=:bottomleft); scatter!(x,y2,label="puntos Interpolados");
savefig(plt1,@sprintf("pInterpolacion1-%i.png",n))
savefig(plt2,@sprintf("pInterpolacion2-%i.png",n))

n=8
x=sort(2*rand(n).-1)
y1=(25*x.^2 .+ 1).^-1
y2=abs.(x)-x/2 - x.^2

pI1=interpolarPolinomios(z,x,y1,m,n)
pI2=interpolarPolinomios(z,x,y2,m,n)
pL1=polinomioLagrange(z,x,y1,m,n)
pL2=polinomioLagrange(z,x,y2,m,n)
pN1=polinomioNewton(z,x,y1,m,n)
pN2=polinomioNewton(z,x,y2,m,n)

plt1=plot(z,[f1,pI1,pL1,pN1],label = [L"\frac{1}{1+25x^2}" "polinomio interpolado" "polinomio de Lagrange" "polinomio de Newton"],legend=:bottomleft); scatter!(x,y1,label="puntos Interpolados");
plt2=plot(z,[f2,pI2,pL2,pN2],label = [L"|x| - \frac{x}{2} - x^2" "polinomio interpolado" "polinomio de Lagrange" "polinomio de Newton"],legend=:topleft); scatter!(x,y2,label="puntos Interpolados");
savefig(plt1,@sprintf("pInterpolacion1-%i.png",n))
savefig(plt2,@sprintf("pInterpolacion2-%i.png",n))




Random.seed!(2022);

function g1(x::AbstractFloat)
    g=[1,exp(-x^2),cospi(x)]
    return g
end

function g2(x::AbstractFloat)
    g=[x,exp(-x^2),cospi(2*x)]
    return g
end

n=3
m=5
z=sort(2*rand(m).-1)
y1=(25*z.^2 .+ 1).^-1
y2=abs.(z)-z/2 - z.^2

pR1=minimosCuadradosPolinomios(z,y1,m,1)
pP1=minimosCuadradosPolinomios(z,y1,m,6)
pG1=minimosCuadrados(g1,z,y1,m,n)
pR2=minimosCuadradosPolinomios(z,y2,m,1)
pP2=minimosCuadradosPolinomios(z,y2,m,6)
pG2=minimosCuadrados(g2,z,y2,m,n)

plt1=plot(z,[y1,pR1,pP1,pG1],label = [L"\frac{1}{1+25x^2}" L"\alpha_0 + \alpha_1 x" L"p_6(x)" L"\alpha^{\top}g_1(x)"]);
plt2=plot(z,[y2,pR2,pP2,pG2],label = [L"|x| - \frac{x}{2} - x^2" L"\alpha_0 + \alpha_1 x" L"p_6(x)" L"\alpha^{\top}g_2(x)"]);
savefig(plt1,@sprintf("pMinimosCuadrados1-%i.png",m))
savefig(plt2,@sprintf("pMinimosCuadrados2-%i.png",m))


m=100
z=sort(2*rand(m).-1)
y1=(25*z.^2 .+ 1).^-1
y2=abs.(z)-z/2 - z.^2

pR1=minimosCuadradosPolinomios(z,y1,m,1)
pP1=minimosCuadradosPolinomios(z,y1,m,6)
pG1=minimosCuadrados(g1,z,y1,m,n)
pR2=minimosCuadradosPolinomios(z,y2,m,1)
pP2=minimosCuadradosPolinomios(z,y2,m,6)
pG2=minimosCuadrados(g2,z,y2,m,n)

plt1=plot(z,[y1,pR1,pP1,pG1],label = [L"\frac{1}{1+25x^2}" L"\alpha_0 + \alpha_1 x" L"p_6(x)" L"\alpha^{\top}g_1(x)"]);
plt2=plot(z,[y2,pR2,pP2,pG2],label = [L"|x| - \frac{x}{2} - x^2" L"\alpha_0 + \alpha_1 x" L"p_6(x)" L"\alpha^{\top}g_2(x)"]);
savefig(plt1,@sprintf("pMinimosCuadrados1-%i.png",m))
savefig(plt2,@sprintf("pMinimosCuadrados2-%i.png",m))


m=1000
z=sort(2*rand(m).-1)
y1=(25*z.^2 .+ 1).^-1
y2=abs.(z)-z/2 - z.^2

pR1=minimosCuadradosPolinomios(z,y1,m,1)
pP1=minimosCuadradosPolinomios(z,y1,m,6)
pG1=minimosCuadrados(g1,z,y1,m,n)
pR2=minimosCuadradosPolinomios(z,y2,m,1)
pP2=minimosCuadradosPolinomios(z,y2,m,6)
pG2=minimosCuadrados(g2,z,y2,m,n)

plt1=plot(z,[y1,pR1,pP1,pG1],label = [L"\frac{1}{1+25x^2}" L"\alpha_0 + \alpha_1 x" L"p_6(x)" L"\alpha^{\top}g_1(x)"]);
plt2=plot(z,[y2,pR2,pP2,pG2],label = [L"|x| - \frac{x}{2} - x^2" L"\alpha_0 + \alpha_1 x" L"p_6(x)" L"\alpha^{\top}g_2(x)"]);
savefig(plt1,@sprintf("pMinimosCuadrados1-%i.png",m))
savefig(plt2,@sprintf("pMinimosCuadrados2-%i.png",m))