include("factQR.jl")
using Printf


#Este progrma resuelve el problema de la tarea, solo cambiamos el número de nodos en la siguiente linea si deseamos. Para n=100 toma aproximadamente 2 segundos, para n=1000 toma muchas horas
A,m,n=leerMatriz("M_sys_125x125.txt")
b=leerVector("V_sys_125x1.txt")

@time x=resuelveQR(A,b,m)

