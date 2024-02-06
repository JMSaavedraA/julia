Para correr simplemente ejecutar

julia main.jl

desde la linea de comando en la carpeta contenedora, o desde Visual Studio Code, igualmente con la carpeta contenedora

Se incluye Solvers.jl y problemaEliptico.jl de tareas pasadas pues algunas funciones de ahí que son reutilizadas.

Las funciones nuevas están contenidas en metodoJacobi.jl, que incluye el método de Jacobi para obtener los eigenpares de una matriz, pero main.jl hace uso de ellas para resolver el problema de la tarea. Se puede cambiar el número de nodos para las diferencias finitas directo en el main.

Además se generan las matrices P y lambda de cada uno de los métodos en un archivo .txt