#!/bin/bash

# Compilar los programas
mpicxx ADA05_SendRecv.cpp -o ADA05_SendRecvExec
mpicxx ADA05_Reduce.cpp -o ADA05_ReduceExec
mpicxx ADA05_Allreduce.cpp -o ADA05_AllreduceExec
mpicxx ADA05_Gather.cpp -o ADA05_GatherExec
mpicxx ADA05_Allgather.cpp -o ADA05_AllgatherExec

# Mensaje
echo "Método de comunicación: SendRecv"

# Bucle para ejecutar cada programa con diferentes números de procesos
for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./ADA05_SendRecvExec
done

echo "Método de comunicación: Reduce"

for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./ADA05_ReduceExec
done

echo "Método de comunicación: Allreduce"

for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./ADA05_AllreduceExec
done

echo "Método de comunicación: Gather"

for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./ADA05_GatherExec
done

echo "Método de comunicación: Allgather"

for ((i = 1; i <= 18; i++)); do
    echo "Ejecutando con $i procesos..."
    mpiexec -n $i ./ADA05_AllgatherExec
done