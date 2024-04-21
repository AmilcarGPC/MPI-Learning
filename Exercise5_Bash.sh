#!/bin/bash

# Compilar los programas
mpicxx Exercise5_MPI.cpp -o Exercise5_MPI

# Mensaje
echo "Ejecutando en C++"

# Bucle para ejecutar cada programa con diferentes números de procesos
for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./Exercise5_MPI
done