include("Solvers.jl") #Reutilizamos métodos de la Tarea 03
include("Cholesky.jl") #Estos son los métodos construidos para esta Tarea
using Printf

#El número de nodos del problema
nodos=100
println(@sprintf("Problema Elíptico Unidimensional con %i nodos",nodos))

#La solución exacta es z
xs=Base._linspace(0.0,1.0,nodos)
z=xs.^2+xs

#Metodo Analítico (Factorización de Cholesky)
println("Solución por Método Analítico con Factorización de Cholesky")
println("El tiempo de máquina es:")
@time x1=problemaElipticoDiagonal(nodos)
eAbsoluto1=norma2(x1-z,nodos)
eRelativo1=norma2(x1-z,nodos)/norma2(z,nodos)
println(@sprintf("El error absoluto es %.15e",eAbsoluto1))
println(@sprintf("El error relativo es %.15e",eRelativo1))

#Metodo de Gauss-Seidel
println("Solución por Método de Gauss-Seidel")
println("El tiempo de máquina es:")
@time x2,Cercania2,k2=problemaElipticoGaussSeidel(nodos)
eAbsoluto2=norma2(x2-z,nodos)
eRelativo2=norma2(x2-z,nodos)/norma2(z,nodos)
println(@sprintf("El error absoluto es %.15e",eAbsoluto2))
println(@sprintf("El error relativo es %.15e",eRelativo2))
println(@sprintf("Se realizaron  %i iteraciones",k2))
println(@sprintf("Criterio de paro calculado %.15e",Cercania2[k2]))

#Metodo de Jacobi
println("Solución por Método de Jacobi")
println("El tiempo de máquina es:")
@time x3,Cercania3,k3=problemaElipticoJacobi(nodos)
eAbsoluto3=norma2(x3-z,nodos)
eRelativo3=norma2(x3-z,nodos)/norma2(z,nodos)
println(@sprintf("El error absoluto es %.15e",eAbsoluto3))
println(@sprintf("El error relativo es %.15e",eRelativo3))
println(@sprintf("Se realizaron  %i iteraciones",k3))
println(@sprintf("Criterio de paro calculado %.15e",Cercania3[k3]))