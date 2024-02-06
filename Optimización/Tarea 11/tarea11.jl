#Jose Miguel Saavedra Aguilar
#CIMAT Master in Applied Mathematics
#Optimization. Homework 11


include("funcionAjuste.jl")
include("Levenberg-Marquart.jl")

println("Primer punto inicial\n")

z0 = [15,0.6,0.0];

R = residual(z0);

fz0 = 1/2 * R'*R;

z1, f, k, nP, sol = LevenbergMarquart(residual,jacobianaResidual,z0,0.001);

showFitResults(z1,f,k,nP,sol)

plt1 = scatter(x,y,label="puntos originales"); plot!(x,y+residual(z1),width = 3);

savefig(plt1,"grafica1.png");

println("\nSegundo punto inicial\n")

z0 = [15,1.0,0.0];

z2, f, k, nP, sol = LevenbergMarquart(residual,jacobianaResidual,z0,100);

showFitResults(z2,f,k,nP,sol)

plt2 = scatter(x,y,label="puntos originales"); plot!(x,y+residual(z2),width = 3);

savefig(plt2,"grafica2.png");

println("\nTercer punto inicial\n")

z0 = [15,0.6,1.6];

z3, f, k, nP, sol = LevenbergMarquart(residual,jacobianaResidual,z0,1e+1);

showFitResults(z3,f,k,nP,sol)

plt3 = scatter(x,y,label="puntos originales"); plot!(x,y+residual(z3),width = 3);

savefig(plt3,"grafica3.png");