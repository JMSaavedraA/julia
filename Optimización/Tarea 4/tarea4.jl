#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics
#Optimization. Homework 4

using Plots, LaTeXStrings, DelimitedFiles
include("funcionesPrueba.jl")
include("optimizadores.jl")

Random.seed!(2023)
# Ejercicio 1
# n=2
#Start with given x_0 for both steepest descent and Newton's point directions.
x0=[-1.2 ; 1]
xPF1,gPF1,FPF1,GPF1=pasoFijo(Rosenbrock2,gradienteRosenbrock2,x0,50000,1e-8,0.0019)
xN1,gN1,FN1,GN1=descensoNewton(Rosenbrock2,gradienteRosenbrock2,hessianaRosenbrock2,x0,5000,1e-8)

#Continue with a random starting point for both steepest descent and Newton's point directions.
x0=rand(2)
xPF2,gPF2,FPF2,GPF2=pasoFijo(Rosenbrock2,gradienteRosenbrock2,x0,50000,1e-8,0.0019)
xN2,gN2,FN2,GN2=descensoNewton(Rosenbrock2,gradienteRosenbrock2,hessianaRosenbrock2,x0,5000,1e-8)

#Save a plot of (k,f(x_k)) and one of (k,∇f(x_k))
plt=plot(FPF1,size=[800,600],label=L"x_0",lw=3);plot!(FPF2,size=[800,600],label="Aleatorio",lw=3)
savefig(plt,"FPF1.png")
plt=plot(FN1,size=[800,600],label=L"x_0",lw=3);plot!(FN2,size=[800,600],label="Aleatorio",lw=3)
savefig(plt,"FN1.png")
plt=plot(GPF1,size=[800,600],label=L"x_0",lw=3);plot!(GPF2,size=[800,600],label="Aleatorio",lw=3)
savefig(plt,"GPF1.png")
plt=plot(GN1,size=[800,600],label=L"x_0",lw=3);plot!(GN2,size=[800,600],label="Aleatorio",lw=3)
savefig(plt,"GN1.png")

# n=100
#Start with given x_0 for both steepest descent and Newton's point directions.
x0=ones(100)
x0[1]=-1.2
x0[99]=-1.2
xPF3,gPF3,FPF3,GPF3=pasoFijo(Rosenbrock100,gradienteRosenbrock100,x0,50000,1e-8,0.0011)
xN3,gN3,FN3,GN3=descensoNewton(Rosenbrock100,gradienteRosenbrock100,hessianaRosenbrock100,x0,5000,1e-8)

#Continue with a random starting point for both steepest descent and Newton's point directions.
x0=rand(100)
xPF4,gPF4,FPF4,GPF4=pasoFijo(Rosenbrock100,gradienteRosenbrock100,x0,50000,1e-8,0.0011)
xN4,gN4,FN4,GN4=descensoNewton(Rosenbrock100,gradienteRosenbrock100,hessianaRosenbrock100,x0,5000,1e-8)

#Save a plot of (k,f(x_k)) and one of (k,∇f(x_k))
plt=plot(FPF3,size=[800,600],label=L"x_0",lw=3);plot!(FPF4,size=[800,600],label="Aleatorio",lw=3)
savefig(plt,"FPF3.png")
plt=plot(FN3,size=[800,600],label=L"x_0",lw=3);plot!(FN4,size=[800,600],label="Aleatorio",lw=3)
savefig(plt,"FN3.png")
plt=plot(GPF3,size=[800,600],label=L"x_0",lw=3);plot!(GPF4,size=[800,600],label="Aleatorio",lw=3)
savefig(plt,"GPF3.png")
plt=plot(GN3,size=[800,600],label=L"x_0",lw=3);plot!(GN4,size=[800,600],label="Aleatorio",lw=3)
savefig(plt,"GN3.png")

#Ejercicio 2
#Start with given x_0 for both steepest descent and Newton's point directions.
x0=[-3,-1,-3,-1]
xPF5,gPF5,FPF5,GPF5=pasoFijo(Wood,gradienteWood,x0,50000,1e-8,0.000421)
xN5,gN5,FN5,GN5=descensoNewton(Wood,gradienteWood,hessianaWood,x0,5000,1e-8)

#Continue with a random starting point for both steepest descent and Newton's point directions.
x0=rand(4)
xPF6,gPF6,FPF6,GPF6=pasoFijo(Wood,gradienteWood,x0,50000,1e-8,0.001)
xN6,gN6,FN6,GN6=descensoNewton(Wood,gradienteWood,hessianaWood,x0,5000,1e-8)

#Save a plot of (k,f(x_k)) and one of (k,∇f(x_k))
plt=plot(FPF5,size=[800,600],label=L"x_0",lw=3);plot!(FPF6,size=[800,600],label="Aleatorio",lw=3)
savefig(plt,"FPF5.png")
plt=plot(FN5,size=[800,600],label=L"x_0",lw=3);plot!(FN6,size=[800,600],label="Aleatorio",lw=3)
savefig(plt,"FN5.png")
plt=plot(GPF5,size=[800,600],label=L"x_0",lw=3);plot!(GPF6,size=[800,600],label="Aleatorio",lw=3)
savefig(plt,"GPF5.png")
plt=plot(GN5,size=[800,600],label=L"x_0",lw=3);plot!(GN6,size=[800,600],label="Aleatorio",lw=3)
savefig(plt,"GN5.png")

#Ejercicio 3
global η=randn(128) #Simulate n=128 normally distributed numbers μ=0 and σ=1
global Y=((2*(1:128) .- 129) ./ 127).^2 + η #Calculate y_i. Note this must be done before calling "funImagen" function. funImagen and it's gradient will automatically take it as it's a global variable
global λ::Int #λ can only be an integer value.

# λ=1
λ=1
#Start from x_0=Y for both steepest descent and Newton's point directions.
x0 = copy(Y)
xPF7,gPF7,FPF7,GPF7=pasoFijo(funImagen,gradImagen,x0,50000,1e-10,0.1)
xN7,gN7,FN7,GN7=descensoNewton(funImagen,gradImagen,hessianaImagen,x0,5000,1e-10)

#Continue with a random starting point for both steepest descent and Newton's point directions.
x0=rand(128)
xPF8,gPF8,FPF8,GPF8=pasoFijo(funImagen,gradImagen,x0,5000,1e-10,0.1)
xN8,gN8,FN8,GN8=descensoNewton(funImagen,gradImagen,hessianaImagen,x0,5000,1e-10)


# λ=10
λ=10
#Start from x_0=Y for both steepest descent and Newton's point directions.
x0 = copy(Y)
xPF9,gPF9,FPF9,GPF9=pasoFijo(funImagen,gradImagen,x0,50000,1e-10,0.0195)
xN9,gN9,FN9,GN9=descensoNewton(funImagen,gradImagen,hessianaImagen,x0,5000,1e-10)

#Continue with a random starting point for both steepest descent and Newton's point directions.
x0=rand(128)
xPF10,gPF10,FPF10,GPF10=pasoFijo(funImagen,gradImagen,x0,50000,1e-10,0.0195)
xN10,gN10,FN10,GN10=descensoNewton(funImagen,gradImagen,hessianaImagen,x0,5000,1e-10)


# λ=1000
λ=1000
#Start from x_0=Y for both steepest descent and Newton's point directions.
x0 = copy(Y)
xPF11,gPF11,FPF11,GPF11=pasoFijo(funImagen,gradImagen,x0,50000,1e-10,0.0002)
xN11,gN11,FN11,GN11=descensoNewton(funImagen,gradImagen,hessianaImagen,x0,5000,1e-10)

#Continue with a random starting point for both steepest descent and Newton's point directions.
x0=rand(128)
xPF12,gPF12,FPF12,GPF12=pasoFijo(funImagen,gradImagen,x0,50000,1e-10,0.0002)
xN12,gN12,FN12,GN12=descensoNewton(funImagen,gradImagen,hessianaImagen,x0,5000,1e-10)

#Scatter plot t_i, y_i, x_i and save to working folder
sct=scatter([((2*(1:128) .- 129) ./ 127).^2,Y,xPF7,xPF9,xPF11], label=[L"t_i" L"y_i" L"\lambda = 1" L"\lambda = 10" L"\lambda = 1000"],color_palette=:seaborn_bright6)
savefig(sct,"sctPF1x0.png")
sct=scatter([((2*(1:128) .- 129) ./ 127).^2,Y,xN7,xN9,xN11], label=[L"t_i" L"y_i" L"\lambda = 1" L"\lambda = 10" L"\lambda = 1000"],color_palette=:seaborn_bright6)
savefig(sct,"sctN1x0.png")

sct=scatter([((2*(1:128) .- 129) ./ 127).^2,Y,xPF8,xPF10,xPF12], label=[L"t_i" L"y_i" L"\lambda = 1" L"\lambda = 10" L"\lambda = 1000"],color_palette=:seaborn_bright6)
savefig(sct,"sctPF1Random.png")
sct=scatter([((2*(1:128) .- 129) ./ 127).^2,Y,xN8,xN10,xN12], label=[L"t_i" L"y_i" L"\lambda = 1" L"\lambda = 10" L"\lambda = 1000"],color_palette=:seaborn_bright6)
savefig(sct,"sctN1Random.png")

#Now, the case for given data file y.txt
global Y=vec(readdlm("y.txt",',')) #y is read from the file "y.txt". funImagen and it's gradient will automatically take it as it's a global variable

# λ=1
λ=1
#Start from x_0=Y for both steepest descent and Newton's point directions.
x0 = copy(Y)
xPF13,gPF13,FPF13,GPF13=pasoFijo(funImagen,gradImagen,x0,50000,1e-10,0.1)
xN13,gN13,FN13,GN13=descensoNewton(funImagen,gradImagen,hessianaImagen,x0,5000,1e-10)

#Continue with a random starting point for both steepest descent and Newton's point directions.
x0=rand(256)
xPF14,gPF14,FPF14,GPF14=pasoFijo(funImagen,gradImagen,x0,50000,1e-10,0.1)
xN14,gN14,FN14,GN14=descensoNewton(funImagen,gradImagen,hessianaImagen,x0,5000,1e-10)


# λ=10
λ=10
#Start from x_0=Y for both steepest descent and Newton's point directions.
x0 = copy(Y)
xPF15,gPF15,FPF15,GPF15=pasoFijo(funImagen,gradImagen,x0,50000,1e-10,0.0195)
xN15,gN15,FN15,GN15=descensoNewton(funImagen,gradImagen,hessianaImagen,x0,5000,1e-10)

#Continue with a random starting point for both steepest descent and Newton's point directions.
x0=rand(256)
xPF16,gPF16,FPF16,GPF16=pasoFijo(funImagen,gradImagen,x0,50000,1e-10,0.0195)
xN16,gN16,FN16,GN16=descensoNewton(funImagen,gradImagen,hessianaImagen,x0,5000,1e-10)


# λ=1000
λ=1000
#Start from x_0=Y for both steepest descent and Newton's point directions.
x0 = copy(Y)
xPF17,gPF17,FPF17,GPF17=pasoFijo(funImagen,gradImagen,x0,50000,1e-8,0.0002)
xN17,gN17,FN17,GN17=descensoNewton(funImagen,gradImagen,hessianaImagen,x0,50000,1e-8)

#Continue with a random starting point for both steepest descent and Newton's point directions.
x0=rand(256)
xPF18,gPF18,FPF18,GPF18=pasoFijo(funImagen,gradImagen,x0,50000,1e-8,0.0002)
xN18,gN18,FN18,GN18=descensoNewton(funImagen,gradImagen,hessianaImagen,x0,50000,1e-8)

#Scatter plot t_i, y_i, x_i and save to working folder
sct=scatter([Y,xPF13,xPF15,xPF17], label=[L"y_i" L"\lambda = 1" L"\lambda = 10" L"\lambda = 1000"],color_palette=:seaborn_bright6)
savefig(sct,"sctPF2x0.png")
sct=scatter([Y,xN13,xN15,xN17], label=[L"y_i" L"\lambda = 1" L"\lambda = 10" L"\lambda = 1000"],color_palette=:seaborn_bright6)
savefig(sct,"sctN2x0.png")

sct=scatter([Y,xPF14,xPF16,xPF18], label=[L"y_i" L"\lambda = 1" L"\lambda = 10" L"\lambda = 1000"],color_palette=:seaborn_bright6)
savefig(sct,"sctPF2Random.png")
sct=scatter([Y,xN14,xN16,xN18], label=[L"y_i" L"\lambda = 1" L"\lambda = 10" L"\lambda = 1000"],color_palette=:seaborn_bright6)
savefig(sct,"sctN2Random.png")