include("metodosTridiagonales.jl")
include("problemaExamen.jl")
include("Solvers.jl")
include("iteracionsubespacios.jl")
include("Potencia.jl")

m=5000;
n=10
A=crearMatrizExamenPentadiagonal(m);
Aaux=reconstruyeNdiagonal(A,m,2)
#x0,_=factorizaQR(Aaux,m,n)
x0=2*rand(m,n).-1
x0,_=factorizaQR(x0,m,n)
@time PhiPotencia,lambdaPotencia=metodoPotenciaPentadiagonalQR(A,x0,m,1e-8,190000,n)
guardarMatriz(PhiPotencia,m,n,"PhiPotencia.txt")
guardarVector(lambdaPotencia,n,"lambdaPotencia.txt")
@time PhiIteracion,lambdaIteracion=iteracionSubespaciosPotencia(Aaux,m,n,100,1e-8,10000,PhiPotencia)
guardarMatriz(PhiIteracion,m,n,"PhiIteracion.txt")
guardarVector(lambdaIteracion,n,"lambdaIteracion.txt")


L=CholeskyLDLPentadiagonal(A,m)
L=CholeskyLDLPentadiagonal(A,m);
d=zeros(Float64,m);
for i=1:m
    d[i]=L[i,3];
    L[i,3]=1.0;
end
@time PhiInversa,lambdaInversa=metodoPotenciaInversaPentadiagonalQR(L,d,x0,m,1e-8,170000,n)
guardarMatriz(PhiInversa,m,n,"PhiInversa.txt")
guardarVector(lambdaInversa,n,"lambdaInversa.txt")
@time PhiIteracionInversa,lambdaIteracionInversa=iteracionSubespaciosPotenciaInversa(Aaux,m,n,10000,1e-8,10000,PhiInversa)
guardarMatriz(PhiIteracionInversa,m,n,"PhiIteracionInversa.txt")
guardarVector(lambdaIteracionInversa,n,"lambdaIteracionInversa.txt")