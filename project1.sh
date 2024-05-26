#!/bin/bash
OUTPUT_FILE="projectOut.txt"

function encontrar_combinaciones {
    x=$1
    for ((i = 1; i <= x; i++)); do
        for ((j = 1; j <= x; j++)); do
            if [ $((i * j)) -eq $x ]; then
                mpiexec -n $x ./run_test $i $j >> "$OUTPUT_FILE"
            fi
        done
    done
    echo
}

for ((num = 0; num <= 20; num++)); do
    echo "Número de procesos $num:"
    encontrar_combinaciones $num
done