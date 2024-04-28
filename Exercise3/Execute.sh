#!/bin/bash

# Compilar los programas
mpicxx SendRecv_MPI.cpp -o run_SendRecv_MPI
mpicxx Reduce_MPI.cpp -o run_Reduce_MPI
mpicxx Allreduce_MPI.cpp -o run_Allreduce_MPI
mpicxx Gather_MPI.cpp -o run_Gather_MPI
mpicxx Allgather_MPI.cpp -o run_Allgather_MPI

# Mensaje
echo "Método de comunicación: SendRecv"

# Bucle para ejecutar cada programa con diferentes números de procesos
for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./run_SendRecv_MPI
done

echo "Método de comunicación: Reduce"

for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./run_Reduce_MPI
done

echo "Método de comunicación: Allreduce"

for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./run_Allreduce_MPI
done

echo "Método de comunicación: Gather"

for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./run_Gather_MPI
done

echo "Método de comunicación: Allgather"

for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./run_Allgather_MPI
done