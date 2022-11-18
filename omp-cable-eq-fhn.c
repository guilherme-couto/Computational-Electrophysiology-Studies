/* Plot gnuplot
*   gnuplot plot 'v.dat' w l
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int main(int argc, char *argv[])
{

    int num_threads;

    if (argc != 2)
    {
        fprintf(stderr, "Usage: %s <number of threads>\n", argv[0]);
        exit(1);
    }

    num_threads = atoi(argv[1]);

    if (num_threads <= 0)
    {
        fprintf(stderr, "Number of threads must greater than 0\n");
        exit(1);
    }
    printf("The number of threads is %d\n", num_threads);

    // Infos
    float k, A, alpha, L, epsilon, gamma;
    float d_x, d_t;

    k = 2.0;
    A = 1.0;
    alpha = 0.1;
    L = 150;
    epsilon = 0.005;
    gamma = 2.0;

    d_x = 1.0;
    d_t = 0.1;

    int T, M, N;

    T = 80000;

    M = T / d_t; // number of points in time
    N = L / d_x; // number of points in space

    printf("Number of points in time = %d\n", M);
    printf("Number of points in space = %d\n", N);

    // Matrix for v and w
    double *V;
    double *W;
    V = malloc(M * N * sizeof(double));
    W = malloc(M * N * sizeof(double));

    // Initial conditions
    for (int i = 0; i < M * N; i++)
    {
        if (i < N / 10)
        {
            V[i] = 0.3;
        }
        else
        {
            V[i] = 0.0;
        }
        W[i] = 0.0;
    }

    printf("\n");

    // Check time
    double start, finish, elapsed;
    start = omp_get_wtime();

    // Explicit Euler
    for (int n = 0; n < M - 1; n++)
    {
        int i;

        # pragma omp parallel for num_threads(num_threads) default(none) \
        private(i) shared(n, V, W, N, M, k, A, alpha, epsilon, gamma, d_x, d_t)
        
        for (i = 1; i < N - 1; i++)
        {
            V[(n + 1) * N + i] = V[n * N + i] + d_t * ((k * (V[n * N + (i + 1)] - 2 * V[n * N + i] + V[n * N + (i - 1)]) / d_x / d_x) + (A * V[n * N + i] * (1 - V[n * N + i]) * (V[n * N + i] - alpha)) - W[n * N + i]);
            W[(n + 1) * N + i] = W[n * N + i] + d_t * epsilon * (V[(n + 1) * N + i] - gamma * W[n * N + i]);
        }

        // Condition
        V[(n + 1) * N + 0] = V[(n + 1) * N + 1];
        V[(n + 1) * N + (N - 1)] = V[(n + 1) * N + (N - 2)];
    }

    finish = omp_get_wtime();
    elapsed = finish - start;
    printf("Elapsed time = %e seconds\n", elapsed);

    // Print
    for (int i = M - 4; i < M; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            printf("%f ", V[i * N + j]);
        }
        printf("\n");
    }

    FILE *fp = NULL;
    fp = fopen("v.dat", "w");
    for (int i = 0; i < N; i++)
    {
        double a = i;
        fprintf(fp, "%lf\t %lf\n", a, V[1000*N+i]);
    }

    free(V);
    free(W);

    return 0;
}