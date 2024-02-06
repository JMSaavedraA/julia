include("Solvers.jl")
include("rayleigh.jl")
using Printf

#Este progrma resuelve el problema 3 de la tarea, se puede cambiar la matriz en la siguiente linea
A,m,n=leerMatriz("Eigen_50x50.txt")

#Iniciamos en un vector de ceros y uno en la entrada k. Así mismo, lambda aproximada es 10*k. Cambiar k para aproximar otros eigenpares.
V=zeros(Float64,m,n)
Lambda=zeros(Float64,m)
@time for k=1:m
    vAprox=zeros(Float64,m)
    vAprox[k]=1
    lambdaAprox=3.0*k
    V[:,k],Lambda[k]=eigenRayleigh(A,vAprox,lambdaAprox,m,1e-10,100)
end

guardarMatriz(V,m,n,@sprintf("PhiRayleigh%i.txt",m))
guardarVector(Lambda,n,@sprintf("lambdaRayleigh%i.txt",m))

eAbsoluto=A*V-V.*Lambda';

for i=1:n
    println(@sprintf("El error de aproximación del %iº eigenpar es %g",i,norma2(eAbsoluto[:,i],m)))
end