include("factQR.jl")
using Printf


#Este progrma resuelve el problema 4 de la tarea, solo cambiamos el sistema resuelto en las siguientes lineas.
A,m,n=leerMatriz("M_sys_125x125.txt")
b=leerVector("V_sys_125x1.txt")

@time x=resuelveQR(A,b,m)

guardarVector(x,m,@sprintf("xQR%i.txt",m))

println(@sprintf("El error relativo de aproximación es %g" ,norma2(A*x-b,m)/norma2(x,m)))
