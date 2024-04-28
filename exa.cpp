#include <mpi.h>
#include <iostream>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double** bigMatrix;
    double** smallMatrix;

    if (rank == 1) {
        bigMatrix = new double*[10]; // Matriz grande de 10x10
        for (int i = 0; i < 10; i++) {
            bigMatrix[i] = new double[10];
            for (int j = 0; j < 10; j++) {
                bigMatrix[i][j] = 0.0;
            }
        }
    }

    if (rank == 0) {
        smallMatrix = new double*[3]; // Matriz chica de 3x3
        for (int i = 0; i < 3; i++) {
            smallMatrix[i] = new double[3];
        }

        // Inicializar la matriz chica con los valores deseados
        int value = 1;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                smallMatrix[i][j] = value;
                value++;
            }
        }

        // Enviar la primera fila de la matriz chica al proceso 1
        MPI_Send(smallMatrix[0], 3, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    }

    if (rank == 1) {
        // Recibir la primera fila de la matriz chica del proceso 0
        smallMatrix = new double*[1];
        smallMatrix[0] = new double[3];
        MPI_Recv(smallMatrix[0], 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Incorporar la primera fila de la matriz chica en la matriz grande
        for (int j = 0; j < 3; j++) {
            bigMatrix[0][j] = smallMatrix[0][j];
        }

        // Imprimir la matriz grande
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                std::cout << bigMatrix[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    if (rank == 0) {
        for (int i = 0; i < 3; i++) {
            delete[] smallMatrix[i];
        }
        delete[] smallMatrix;
    }

    if (rank == 1) {
        for (int i = 0; i < 10; i++) {
            delete[] bigMatrix[i];
        }
        delete[] smallMatrix[0];
        delete[] smallMatrix;
        delete[] bigMatrix;
    }

    MPI_Finalize();
    return 0;
}
