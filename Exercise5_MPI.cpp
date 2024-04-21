#include <stdio.h>
#include <math.h>

// Definición de la función
double f(double x, double y) {
    return 1 + 6 * x * y * y;
}

// Función para calcular la integral doble utilizando la regla del trapecio
double integrate(double (*f)(double, double), double x0, double x1, double y0, double y1, int nx, int ny) {
    double hx = (x1 - x0) / nx;
    double hy = (y1 - y0) / ny;
    double integral = 0.0;
    
    // Integración en x y y utilizando la regla del trapecio
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            double x = x0 + i * hx;
            double y = y0 + j * hy;
            double weight = 1.0;
            if (i == 0 || i == nx) weight *= 0.5;
            if (j == 0 || j == ny) weight *= 0.5;
            integral += weight * f(x, y);
        }
    }
    integral *= hx * hy;
    return integral;
}

int main() {
    // Definición de los límites de integración
    double x0 = 0.0, x1 = 2.0;
    double y0 = -1.0, y1 = 1.0;
    // Número de divisiones en x e y
    int nx = 100, ny = 100;

    // Calculando la integral doble
    double result = integrate(f, x0, x1, y0, y1, nx, ny);
    
    printf("El resultado de la integral doble es: %f\n", result);
    
    return 0;
}
