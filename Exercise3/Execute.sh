#!/bin/bash

# Compilar los programas
mpicxx SendRecv_MPI.cpp -o SendRecv_MPIExec
mpicxx Reduce_MPI.cpp -o Reduce_MPIExec
mpicxx Allreduce_MPI.cpp -o Allreduce_MPIExec
mpicxx Gather_MPI.cpp -o Gather_MPIExec
mpicxx Allgather_MPI.cpp -o Allgather_MPIExec

# Mensaje
echo "Método de comunicación: SendRecv"

# Bucle para ejecutar cada programa con diferentes números de procesos
for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./SendRecv_MPIExec
done

echo "Método de comunicación: Reduce"

for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./Reduce_MPIExec
done

echo "Método de comunicación: Allreduce"

for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./Allreduce_MPIExec
done

echo "Método de comunicación: Gather"

for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./Gather_MPIExec
done

echo "Método de comunicación: Allgather"

for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./Allgather_MPIExec
done