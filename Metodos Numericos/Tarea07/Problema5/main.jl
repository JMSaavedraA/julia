include("gradienteConjugado.jl")
using Printf


#Este progrma resuelve el problema 5 de la tarea, solo cambiamos el sistema resuelto en las siguientes lineas.
A,m,n=leerMatriz("M_sys_125x125.txt")
b=leerVector("V_sys_125x1.txt")

@time x=gradienteConjugado(A,b,m,1e-10,3*m)

guardarVector(x,m,@sprintf("xCG%i.txt",m))

println(@sprintf("El error relativo de aproximación es %g" ,norma2(A*x-b,m)/norma2(x,m)))