#Importamos las funciones construidas para la Tarea 2
include("Raices.jl")

#Definimos las funciones de nuestro interés. Para el inciso a, podemos modificar el coeficiente en la siguiente linea para hacer las perturbaciones

coeficiente=0.9

function fun1(x)
    y=coeficiente*x^3-21*x^2+120*x-100;
    return y
end
function fun1Prima(x)
    y=coeficiente*3*x^2-42*x+120;
    return y
end

function fun3(x)
    y=log(x^2+1) - exp(0.4*x)*cospi(x);
    return y
end
function fun3Prima(x)
    y=(2*x)/(x^2+1) - exp(0.4*x)*(0.4*cospi(x)-pi*sinpi(x));
    return y
end

#A continuación, se corren los algoritmos para los diferentes puntos de inicio/intervalos de interés que se mencionan en el PDF. Nótese que no todos los valores del coeficiente tendrán convergencia para la Bisección, como se vio en el PDF

x1=Bisección(fun1,0.0,5.0,100,1e-5);
x2=newtonRaphson(fun1,fun1Prima,0.0,100,1e-6);
x3=newtonSinDerivada(fun1,0.0,100,1e-6,1e-8);
x4=Bisección(fun1,5.0,10.0,100,1e-5);
x5=newtonRaphson(fun1,fun1Prima,9.0,100,1e-6);
x6=newtonSinDerivada(fun1,9.0,100,1e-6,1e-8);
x1=Bisección(fun1,10.0,15.0,100,1e-5);
x8=newtonRaphson(fun1,fun1Prima,11.0,100,1e-6);
x9=newtonSinDerivada(fun1,11.0,100,1e-6,1e-8);

println()

#Ahora se consideran el inciso c los puntos comparados en el PDF
println("Para i=0")
x01=Bisección(fun3,-0.95968843,0.04031157172,100,1e-5);
x02=newtonRaphson(fun3,fun3Prima,-0.5,100,1e-6);
x03=newtonSinDerivada(fun3,-0.5,100,1e-6,1e-8);
println()

println("Para i=1")
x11=Bisección(fun3,0.040311571725783,1.040311571725783,100,1e-5);
x12=newtonRaphson(fun3,fun3Prima,0.5,100,1e-6);
x13=newtonSinDerivada(fun3,0.5,100,1e-6,1e-8);
println()

println("Para i=5")
x51=Bisección(fun3,4.040311571725783,5.040311571725783,100,1e-5);
x52=newtonRaphson(fun3,fun3Prima,4.5,100,1e-6);
x53=newtonSinDerivada(fun3,4.5,100,1e-6,1e-8);
println()

println("Para i=20")
x201=Bisección(fun3,19.040311571725783,20.040311571725783,100,1e-5);
x202=newtonRaphson(fun3,fun3Prima,19.5,100,1e-6);
x203=newtonSinDerivada(fun3,19.5,100,1e-6,1e-8);
println()

println("Para i=50")
x501=Bisección(fun3,49.040311571725783,50.040311571725783,100,1e-5);
x502=newtonRaphson(fun3,fun3Prima,49.5,100,1e-6);
x503=newtonSinDerivada(fun3,49.5,100,1e-6,1e-8);