#!/bin/bash

# Compilar los programas
mpicxx Exercise4_MPI.cpp -o Exercise4_MPI
mpif77 pi_MPI.F90 -o pi_MPI

# Mensaje
echo "Ejecutando en C++"

# Bucle para ejecutar cada programa con diferentes números de procesos
for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./Exercise4_MPI
done

echo "Ejecutando en Fortran"

for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./pi_MPI 
done