#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics
#Optimization. Homework 5

using Plots, LaTeXStrings, DelimitedFiles
include("funcionesPrueba.jl")
include("optimizadores.jl")

Random.seed!(2023)
# Rosenbrock's function
# n=2
#Start with given x_0 for both steepest descent and Newton's point directions.
x0=[-1.2 ; 1]
xT1,gT1,FT1,GT1=backtrackingDescent(Rosenbrock2,gradienteRosenbrock2,x0,50000,1e-8,0.5)
xB1,gB1,FB1,GB1=bisectionDescent(Rosenbrock2,gradienteRosenbrock2,x0,50000,1e-8)
xZ1,gZ1,FZ1,GZ1=zoomDescent(Rosenbrock2,gradienteRosenbrock2,x0,50000,1e-8)

#Continue with a random starting point for both steepest descent and Newton's point directions.
x0=rand(2)
xT2,gT2,FT2,GT2=backtrackingDescent(Rosenbrock2,gradienteRosenbrock2,x0,50000,1e-8,0.5)
xB2,gB2,FB2,GB2=bisectionDescent(Rosenbrock2,gradienteRosenbrock2,x0,50000,1e-8)
xZ2,gZ2,FZ2,GZ2=zoomDescent(Rosenbrock2,gradienteRosenbrock2,x0,50000,1e-8)


#Save a plot of (k,f(x_k)) and one of (k,∇f(x_k))
plt=plot(FT1,size=[800,600],label=L"x_0",lw=3);plot!(FT2,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(FT1),maximum(FT2)]),100])])
savefig(plt,"FBacktracking1.png")
plt=plot(FB1,size=[800,600],label=L"x_0",lw=3);plot!(FB2,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(FB1),maximum(FB2)]),100])])
savefig(plt,"FBiseccion1.png")
plt=plot(FZ1,size=[800,600],label=L"x_0",lw=3);plot!(FZ2,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(FZ1),maximum(FZ2)]),100])])
savefig(plt,"FZoom1.png")
plt=plot(GT1,size=[800,600],label=L"x_0",lw=3);plot!(GT2,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(GT1),maximum(GT2)]),100])])
savefig(plt,"GBacktracking1.png")
plt=plot(GB1,size=[800,600],label=L"x_0",lw=3);plot!(GB2,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(GB1),maximum(GB2)]),100])])
savefig(plt,"GBiseccion1.png")
plt=plot(GZ1,size=[800,600],label=L"x_0",lw=3);plot!(GZ2,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(GZ1),maximum(GZ2)]),100])])
savefig(plt,"GZoom1.png")

# n=100
#Start with given x_0 for both steepest descent and Newton's point directions.
x0=ones(100)
x0[1]=-1.2
x0[99]=-1.2
xT3,gT3,FT3,GT3=backtrackingDescent(Rosenbrock100,gradienteRosenbrock100,x0,50000,1e-8,0.5)
xB3,gB3,FB3,GB3=bisectionDescent(Rosenbrock100,gradienteRosenbrock100,x0,50000,1e-8)
xZ3,gZ3,FZ3,GZ3=zoomDescent(Rosenbrock100,gradienteRosenbrock100,x0,50000,1e-8)

#Continue with a random starting point for both steepest descent and Newton's point directions.
x0=rand(100)
xT4,gT4,FT4,GT4=backtrackingDescent(Rosenbrock100,gradienteRosenbrock100,x0,50000,1e-8,0.5)
xB4,gB4,FB4,GB4=bisectionDescent(Rosenbrock100,gradienteRosenbrock100,x0,50000,1e-8)
xZ4,gZ4,FZ4,GZ4=zoomDescent(Rosenbrock100,gradienteRosenbrock100,x0,50000,1e-8)

#Save a plot of (k,f(x_k)) and one of (k,∇f(x_k))
plt=plot(FT3,size=[800,600],label=L"x_0",lw=3);plot!(FT4,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(FT3),maximum(FT4)]),100])])
savefig(plt,"FBacktracking3.png")
plt=plot(FB3,size=[800,600],label=L"x_0",lw=3);plot!(FB4,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(FB3),maximum(FB4)]),100])])
savefig(plt,"FBiseccion3.png")
plt=plot(FZ3,size=[800,600],label=L"x_0",lw=3);plot!(FZ4,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(FZ3),maximum(FZ4)]),100])])
savefig(plt,"FZoom3.png")
plt=plot(GT3,size=[800,600],label=L"x_0",lw=3);plot!(GT4,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(GT3),maximum(GT4)]),100])])
savefig(plt,"GBacktracking3.png")
plt=plot(GB3,size=[800,600],label=L"x_0",lw=3);plot!(GB4,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(GB3),maximum(GB4)]),100])])
savefig(plt,"GBiseccion3.png")
plt=plot(GZ3,size=[800,600],label=L"x_0",lw=3);plot!(GZ4,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(GZ3),maximum(GZ4)]),100])])
savefig(plt,"GZoom3.png")

# Wood Function
x0=[-3.0,-1.0,-3.0,-1.0]
xT5,gT5,FT5,GT5=backtrackingDescent(Wood,gradienteWood,x0,50000,1e-8,0.5)
xB5,gB5,FB5,GB5=bisectionDescent(Wood,gradienteWood,x0,50000,1e-8)
xZ5,gZ5,FZ5,GZ5=zoomDescent(Wood,gradienteWood,x0,50000,1e-8)

#Continue with a random starting point for both steepest descent and Newton's point directions.
x0=rand(4)
xT6,gT6,FT6,GT6=backtrackingDescent(Wood,gradienteWood,x0,50000,1e-8,0.5)
xB6,gB6,FB6,GB6=bisectionDescent(Wood,gradienteWood,x0,50000,1e-8)
xZ6,gZ6,FZ6,GZ6=zoomDescent(Wood,gradienteWood,x0,50000,1e-8)

#Save a plot of (k,f(x_k)) and one of (k,∇f(x_k))
plt=plot(FT5,size=[800,600],label=L"x_0",lw=3);plot!(FT6,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(FT5),maximum(FT6)]),100])])
savefig(plt,"FBacktracking5.png")
plt=plot(FB5,size=[800,600],label=L"x_0",lw=3);plot!(FB6,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(FB5),maximum(FB6)]),100])])
savefig(plt,"FBiseccion5.png")
plt=plot(FZ5,size=[800,600],label=L"x_0",lw=3);plot!(FZ6,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(FZ5),maximum(FZ6)]),100])])
savefig(plt,"FZoom5.png")
plt=plot(GT5,size=[800,600],label=L"x_0",lw=3);plot!(GT6,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(GT5),maximum(GT6)]),100])])
savefig(plt,"GBacktracking5.png")
plt=plot(GB5,size=[800,600],label=L"x_0",lw=3);plot!(GB6,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(GB5),maximum(GB6)]),100])])
savefig(plt,"GBiseccion5.png")
plt=plot(GZ5,size=[800,600],label=L"x_0",lw=3);plot!(GZ6,size=[800,600],label="Aleatorio",lw=3,ylims=[0,minimum([1.1*maximum([maximum(GZ5),maximum(GZ6)]),100])])
savefig(plt,"GZoom5.png")

#Part 4
global η=randn(128) #Simulate n=128 normally distributed numbers μ=0 and σ=1
global Y=((2*(1:128) .- 129) ./ 127).^2 + η #Calculate y_i. Note this must be done before calling "funImagen" function. funImagen and it's gradient will automatically take it as it's a global variable
global λ::Int #λ can only be an integer value.

# λ=1
λ=1
#Start from x_0=Y for both steepest descent and Newton's point directions.
x0 = copy(Y)
xZ7,gZ7,FZ7,GZ7=zoomDescent(funImagen,gradImagen,x0,50000,1e-8)

# λ=10
λ=10
xZ8,gZ8,FZ8,GZ8=zoomDescent(funImagen,gradImagen,x0,50000,1e-8)

# λ=1000
λ=1000
xZ9,gZ9,FZ9,GZ9=zoomDescent(funImagen,gradImagen,x0,50000,1e-8)

#Scatter plot t_i, y_i, x_i and save to working folder
sct=scatter([((2*(1:128) .- 129) ./ 127).^2,Y,xZ7,xZ8,xZ9], label=[L"t_i" L"y_i" L"\lambda = 1" L"\lambda = 10" L"\lambda = 1000"],color_palette=:seaborn_bright6)
savefig(sct,"e4.png")