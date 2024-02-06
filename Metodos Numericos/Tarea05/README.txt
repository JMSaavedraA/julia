Para correr simplemente ejecutar

julia main.jl

desde la linea de comando en la carpeta contenedora, o desde Visual Studio Code, igualmente con la carpeta contenedora

Se incluye Solvers.jl, metodosTridiagonales.jl y problemaEliptico.jl de tareas pasadas pues algunas funciones de ahí que son reutilizadas.

Las funciones nuevas están contenidas en Potencia.jl, estas son los métodos de potencia, potencia inversa, sus versiones generalizadas y sus respectivas versiones tridiagonales, pero main.jl hace uso de ellas. Se puede cambiar el número de nodos para las diferencias finitas directo en el main.

Además se generan las matrices P y lambda de cada uno de los métodos en un archivo .txt