#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics
#Optimization. Partial Exam 1

include("testsHW.jl")
using Random, Printf

Random.seed!(1996)

kD, kE, tD, tE, CD, CE = testTarea(Rosenbrock2,gradienteRosenbrock2,hessianaRosenbrock2,ones(2),30);

println("Tiempos")
@printf("Rosenbrock, %0.5f,  %0.5f \n",tD,tE)

println("Iteraciones")
@printf("Rosenbrock, %0.5f,  %0.5f \n",kD,kE)

@printf("El método de Dogleg resuelve %2i problemas\n",CD)
@printf("El método exacto resuelve %2i problemas\n",CE)