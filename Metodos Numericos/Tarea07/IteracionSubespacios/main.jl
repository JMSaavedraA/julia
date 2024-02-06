include("Solvers.jl")
include("iteracionsubespacios.jl")
using Printf


#Este progrma resuelve el problema de la tarea, solo cambiamos el número de nodos en la siguiente linea si deseamos. Para n=100 toma aproximadamente 2 segundos, para n=1000 toma muchas horas
A, m, _=leerMatriz("Eigen_50x50.txt")
n=5
maxIter=100000
tol=1e-8
iterJacobi=1000000
iterPotencia=10
x0=zeros(Float64,m,n).+1


@time P,D=iteracionSubespaciosPotenciaInversa(A,m,n,maxIter,tol,iterJacobi,iterPotencia,x0)

guardarMatriz(P,m,n,@sprintf("PJacobi%i.txt",m))
guardarVector(D,n,@sprintf("lambdaJacobi%i.txt",m))
