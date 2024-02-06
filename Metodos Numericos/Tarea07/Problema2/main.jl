include("Solvers.jl")
include("iteracionsubespacios.jl")
using Printf


#Este progrma resuelve el problema 2 de la tarea, se puede cambiar la matriz en la siguiente linea
A, m, _=leerMatriz("Eigen_3x3.txt")
#Este es el número de eigenpares que vamos a obtener, se puede cambiar
n=3
maxIter=100000
tol=1e-8
iterJacobi=1000000
#El primer paso hace el método de la potencia inversa algunas veces, por ejemplo 10. Debe ser mayor a 2 para el $x0$ que consideramos
iterPotencia=10
x0=zeros(Float64,m,n).+1

@time P,D=iteracionSubespaciosPotenciaInversa(A,m,n,maxIter,tol,iterJacobi,iterPotencia,x0)

guardarMatriz(P,m,n,@sprintf("PhiIteracion%i.txt",m))
guardarVector(D,n,@sprintf("lambdaIteracion%i.txt",m))

eAbsoluto=A*P-P.*D'

for i=1:n
    println(@sprintf("El error de aproximación del %iº eigenpar es %g",i,norma2(eAbsoluto[:,i],m)))
end