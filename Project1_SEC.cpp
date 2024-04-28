#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define pi 3.14159265359

int main(int argc, char *argv[]) {
    int i, j, n;
    int Nx, Ny, Nt;
    double kx, ky, xI, xF, yI, yF, tF, tI, Dx, Dy, Dt, t;
    double rx, ry, aC, aN, aS, aE, aW;
    double suma, tiempo;
    double *x, *y;
    double **Uold, **Unew;

    // <----------< INICIALIZACIÓN >---------->
    // Parametros
    kx = 1.0; // Termino difusivo en x
    ky = 1.0; // Termino difusivo en y

    xI = 0.0; // Inicio del dominio en x
    xF = 1.0; // Fin del dominio en x
    yI = 0.0; // Inicio del dominio en y
    yF = 1.0; // Fin del dominio en y

    tI = 0.0; // Tiempo inicial
    tF = 0.2; // Tiempo final

    Nt = 200000; // Numero de pasos en t
    Nx = 400;   // Numero de puntos en x
    Ny = 400;   // Numero de puntos en y

    // Discretización
    Dx = (xF - xI) / (Nx - 1);
    Dy = (yF - yI) / (Ny - 1);
    Dt = (tF - tI) / (Nt - 1);
    rx = kx * (Dt / (Dx * Dx));
    ry = ky * (Dt / (Dy * Dy));

    x = (double *)malloc(Nx * sizeof(double));
    y = (double *)malloc(Ny * sizeof(double));
    Uold = (double **)malloc(Nx * sizeof(double *));
    Unew = (double **)malloc(Nx * sizeof(double *));
    for (i = 0; i < Nx; i++) {
        Uold[i] = (double *)malloc(Ny * sizeof(double));
        Unew[i] = (double *)malloc(Ny * sizeof(double));
    }

    // Malla
    for (i = 0; i < Nx; i++)
        x[i] = xI + i * Dx;
    for (j = 0; j < Ny; j++)
        y[j] = yI + j * Dy;

    // Condiciones Iniciales
    for (i = 0; i < Nx; i++)
        for (j = 0; j < Ny; j++)
            Uold[i][j] = 3.0 * sin(pi * x[i] + pi * y[j]) * sin(pi * x[i] + pi * y[j]);

    // Condiciones de frontera
    for (j = 0; j < Ny; j++) {
        Uold[0][j] = 2.0;    // Oeste
        Uold[Nx - 1][j] = 1.0; // Este
    }
    for (i = 0; i < Nx; i++) {
        Uold[i][0] = 1.0;    // Sur
        Uold[i][Ny - 1] = 3.0; // Norte
    }

    // Actualizacion
    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            Unew[i][j] = Uold[i][j];
        }
    }

    // <----------< PROGRAMA PRINCIPAL >---------->
    // MPI: Loop de calculos
    aC = 1 - 2.0 * (rx + ry);
    aE = rx;
    aW = rx;
    aS = ry;
    aN = ry;

    clock_t inicio = clock();
    for (n = 0; n < Nt; n++) {
        // Nuevos valores
        for (i = 1; i < Nx - 1; i++) {
            for (j = 1; j < Ny - 1; j++) {
                Unew[i][j] = aC * Uold[i][j] +
                             aW * Uold[i - 1][j] +
                             aE * Uold[i + 1][j] +
                             aS * Uold[i][j - 1] +
                             aN * Uold[i][j + 1];
            }
        }

        // Actualización
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                Uold[i][j] = Unew[i][j];
            }
        }
    }
    clock_t final = clock();
    // <----------< VERIFICAR RESULTADOS >---------->
    // Imprimir suma de los elementos para verificar
    suma = 0.0;
    for (i = 0; i < Nx; i++) { 
        for (j = 0; j < Ny; j++) {
            suma += Unew[i][j];
        }
    }
    // Tiempo final
    tiempo = (double)(final - inicio) / CLOCKS_PER_SEC;
    printf("Suma total = %f\n", suma);
    printf("tiempo=%f\n", tiempo);

    // <----------< FINALIZACIÓN >---------->
    // Liberar memoria
    free(x);
    free(y);
    for (i = 0; i < Nx; i++) {
        free(Uold[i]);
        free(Unew[i]);
    }
    free(Uold);
    free(Unew);
    return 0;
}