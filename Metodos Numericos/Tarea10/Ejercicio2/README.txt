Ejercicio 2 de la tarea
Para correr se requieren los paquetes Plots y LaTeXStrings en julia. Se pueden instalar con la orden

julia runFirst.jl

desde la carpeta contenedora. Si se han instalado en las tareas anteriores no es necesario.

Para correr simplemente ejecutar

julia main.jl

desde la linea de comando en la carpeta contenedora.

Se puede cambiar la semilla para la aleatorización de los puntos y_k.

La gráfica se guarda como "elementoFinitoComparacion.png" y cada una de las aproximaciones por elementos finitos como "elementoFinito-lambda.txt" donde lambda es el valor de la penalización, tomamos lambda=0.5, 1.5 y 3.5 pero se puede cambiar en el código.