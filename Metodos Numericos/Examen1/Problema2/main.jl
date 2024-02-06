include("metodosTridiagonales.jl")
include("problemaExamen.jl")
include("Solvers.jl")
using Statistics
using Printf
m=5000;
A=crearMatrizExamenPentadiagonal(m);
AI=zeros(Float64,m,m);
b=zeros(Float64,m);
L=CholeskyLDLPentadiagonal(A,m);
d=copy(b);
for i=1:m
    d[i]=L[i,3];
    L[i,3]=1.0;
end
@time for i=1:m
    b[i]=1.0;
    x=resuelveCholeskyPentadiagonalDada(L,d,b,m);
    b[i]=0.0;
    AI[:,i].=x;
end
Aaux=reconstruyeNdiagonal(A,m,2);
I=Identidad(m);
B=abs.(Aaux*AI-I)
eMax=maximum(B);
eAvg=mean(B)
println(@sprintf("El máximo error de la inversa es de %g",eMax))
println(@sprintf("El error promedio de la inversa es de %g",eAvg))
guardarMatriz(AI,m,m,"InversaExamen.txt")