#!/bin/bash

# Grupo en que ejecutar
#SBATCH --partition=C1Mitad1

# Nombre del trabajo
#SBATCH --job-name=Tarea1OE

#Número máximo de tareas asignadas por nodo.
#SBATCH --ntasks-per-node=12

#Número de procesos a ejecutar, nodos * ntasks-per-node
#SBATCH --ntasks=24

#Archivo de log
#SBATCH --output=logTarea1OE.log

#SBATCH --time=00:10:00

julia instances.jl 21 2500-10 2023 aleatorio2500-21.txt
