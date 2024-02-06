#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics
#Optimization. Homework 8

include("testsHW.jl")
using Random, Printf

Random.seed!(2023)

kRN, kRD, kRB, kRF, tRN, tRD, tRB, tRF = testTarea(Rosenbrock100,gradienteRosenbrock100,hessianaRosenbrock100,ones(100),30);
kWN, kWD, kWB, kWF, tWN, tWD, tWB, tWF = testTarea(Wood,gradienteWood,hessianaWood,ones(4),30);
kBN, kBD, kBB, kBF, tBN, tBD, tBB, tBF = testTarea(Branin,gradienteBranin,hessianaBranin,[π,2.275],30);

println("Tiempos")
@printf("Rosenbrock, %0.5f,  %0.5f,  %0.5f,  %0.5f \n",tRD,tRN,tRB,tRF)
@printf("Wood,       %0.5f,  %0.5f,  %0.5f,  %0.5f \n",tWD,tWN,tWB,tWF)
@printf("Branin,     %0.5f,  %0.5f,  %0.5f,  %0.5f \n",tBD,tBN,tBB,tBF)

println("Iteraciones")
@printf("Rosenbrock, %5.5f,  %5.5f,  %5.5f,    %5.5f \n",kRD,kRN,kRB,kRF)
@printf("Wood,       %5.5f,  %5.5f,  %5.5f,    %5.5f \n",kWD,kWN,kWB,kWF)
@printf("Branin,     %5.5f,  %5.5f,  %5.5f,    %5.5f \n",kBD,kBN,kBB,kBF)