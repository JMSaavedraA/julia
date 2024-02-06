#Primero, importamos nuestro optimizador de funciones

include("Optimizador.jl")

#Lo probaremos en la función 2 de la tarea 2 que no tenía raíces

function fun2(x)
    y=2-log(x)/x;
    return y
end

#Lo corremos y tenemos el mínimo de fun2 en R

xM=optimizadorNewton(fun2,1.0,100,1e-8,1e-6);

println(@sprintf("f(x)=%f es el mínimo de f",fun2(xM)))


#Ahora, importamos las funciones para Ax=b
include("Solvers.jl")

#Leemos la matriz "Ejemplo.txt" y el vector "vector.txt" o "A.txt" y "b.txt". Comentar el que no hagamos y descomentar el otro

A,m,n = leerMatriz("Ejemplo.txt")
b=leerVector("vector.txt")

#A,m,n = leerMatriz("A.txt")
#b=leerVector("b.txt")

#Descomentar el método que querramos utilizar y comentar el resto, para no perder el resultado ni hacer operaciones de más

#@time x=resuelveGaussJordan(A,b,m);

#@time x=resuelveLDU(A,b,m);

#@time x=resuelveCrout(A,b,m);

@time x=resuelveDoolittle(A,b,m);

#A continuación, calculamos la norma infinito del error

print(maximum(abs.(A*x-b),dims=1))

#Finalmente, guardamos el resultado x en el archivo x.txt

guardarVector(x,m)