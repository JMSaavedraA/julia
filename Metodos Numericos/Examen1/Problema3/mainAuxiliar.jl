include("metodosTridiagonales.jl")
include("problemaExamen.jl")
include("Solvers.jl")
include("iteracionsubespacios.jl")
include("Potencia.jl")
include("rayleigh.jl")

m=5000;
n=10
A=crearMatrizExamenPentadiagonal(m);
Aaux=reconstruyeNdiagonal(A,m,2)
PhiPotencia,m,n=leerMatriz("PhiPotencia.txt")
lambdaPotencia=leerVector("lambdaPotencia.txt")
PhiIteracion,m,n=leerMatriz("PhiIteracion.txt")
lambdaIteracion=leerVector("lambdaIteracion.txt")
PhiInversa,m,n=leerMatriz("PhiInversa.txt")
lambdaInversa=leerVector("lambdaInversa.txt")
PhiIteracionInversa,m,n=leerMatriz("PhiIteracionInversa.txt")
lambdaIteracionInversa=leerVector("lambdaIteracionInversa.txt")


ePotencia=zeros(Float64,n)
eIteracion=zeros(Float64,n)
eInversa=zeros(Float64,n)
eIteracionInversa=zeros(Float64,n)

for i=1:n
    ePotencia[i]=norma2(Aaux*PhiPotencia[:,i]-lambdaPotencia[i]*PhiPotencia[:,i],m)
    eIteracion[i]=norma2(Aaux*PhiIteracion[:,i]-lambdaIteracion[i]*PhiIteracion[:,i],m)
    eInversa[i]=norma2(Aaux*PhiInversa[:,i]-lambdaInversa[i]*PhiInversa[:,i],m)
    eIteracionInversa[i]=norma2(Aaux*PhiIteracionInversa[:,i]-lambdaIteracionInversa[i]*PhiIteracionInversa[:,i],m)
end