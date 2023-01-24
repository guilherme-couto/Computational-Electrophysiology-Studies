#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>

// Parameters to reproduce midmyocardial cell
double u_o = 0;
double u_u = 1.61;
double theta_v = 0.3;
double theta_w = 0.13;
double theta_vminus = 0.1;
double theta_o = 0.005;
double tau_v1minus = 80;
double tau_v2minus = 1.4506;
double tau_vplus = 1.4506;
double tau_w1minus = 70;
double tau_w2minus = 8;
double k_wminus = 200;
double u_wminus = 0.016;
double tau_wplus = 280;
double tau_fi = 0.078;
double tau_o1 = 410;
double tau_o2 = 7;
double tau_so1 = 91;
double tau_so2 = 0.8;
double k_so = 2.1;
double u_so = 0.6;
double tau_s1 = 2.7342;
double tau_s2 = 4;
double k_s = 2.0994;
double u_s = 0.9087;
double tau_si = 3.3849;
double tau_winf = 0.01;
double w_infstar = 0.5;

double D = 1.171; // +- 0.0221 cm^2/s  human ventricular diffusion coefficient

// Standar Heaviside function
double H(double x, double y)
{
    if (x > y)
    {
        return 1;
    }
    else if (x < y)
    {
        return 0;
    }
    else
        return 0.5;
}

double v_inf_function(double x, double y)
{
    if (x < y)
    {
        return 1;
    }
    else
        return 0;
}

/*
solves Ax = v where A is a tridiagonal matrix consisting of vectors a, b, c
x - initially contains the input vector v, and returns the solution x. indexed from 0 to X - 1 inclusive
X - number of equations (length of vector x)
a - subdiagonal (means it is the diagonal below the main diagonal), indexed from 1 to X - 1 inclusive
b - the main diagonal, indexed from 0 to X - 1 inclusive
c - superdiagonal (means it is the diagonal above the main diagonal), indexed from 0 to X - 2 inclusive

Note 2: We don't check for diagonal dominance, etc.; this is not guaranteed stable
*/
double* solve_tridiagonal_in_place_reusable(double *d, int N, double *a, double *b, double *c)
{
    /* Allocate scratch space. */
    double *cprime = malloc(N * sizeof(double));
    double *dprime = malloc(N * sizeof(double));

    // Variable vector
    double *y = malloc(N * sizeof(double));

    if (!cprime)
    {
        printf("Error: malloc failed.\n");
    }
    cprime[0] = c[0] / b[0];
    dprime[0] = d[0] / b[0];

    /* loop from 1 to N - 1 inclusive */
    for (int ix = 1; ix < N; ix++)
    {
        double m = 1.0 / (b[ix] - a[ix] * cprime[ix - 1]);
        cprime[ix] = c[ix] * m;
        dprime[ix] = (d[ix] - a[ix] * dprime[ix - 1]) * m;
    }

    y[N - 1] = dprime[N - 1];
    /* loop from N - 2 to 0 inclusive, safely testing loop end condition */
    for (int ix = N - 2; ix >= 0; ix--)
        y[ix] = dprime[ix] - (cprime[ix] * y[ix + 1]);

    /* free scratch space */
    free(cprime);
    free(dprime);

    return y;
}

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

    // Discretization
    int T = 600;
    double L = 300;
    double delta_x = 1;
    double delta_t = 0.001;
    double delta_t_exp = 0.001;

    int M = T / delta_t;               // number of points in time
    int N = L / delta_x;               // number of points in space
    int M_exp = delta_t / delta_t_exp; // number of points in space for explicit method

    printf("Number of points in time = %d\n", M);
    printf("Number of points in space = %d\n", N);
    printf("Number of points in time for explicit method = %d\n", M_exp);

    // Parameters
    double *U, *V, *W, *S;
    U = malloc(M * N * sizeof(double));
    V = malloc(N * sizeof(double));
    W = malloc(N * sizeof(double));
    S = malloc(N * sizeof(double));

    double J_fi, J_so, J_si, J, I_app;
    double tau_vminus, tau_wminus, tau_so, tau_s, tau_o;
    double v_inf, w_inf;
    double du_dt, dv_dt, ds_dt, dw_dt;

    // Initial conditions
    for (int i = 0; i < N; i++)
    {
        U[i] = 0;
        V[i] = 1;
        W[i] = 1;
        S[i] = 0;
    }

    // Subdiagonal, diagonal and superdiagonal of the tridiagonal matrix (Implicit Method)
    double alpha = (D * delta_t) / (delta_x * delta_x);
    double *a, *b, *c, *d, *x;
    a = malloc(N * sizeof(double)); // subdiagonal
    b = malloc(N * sizeof(double)); // diagonal
    c = malloc(N * sizeof(double)); // superdiagonal
    d = malloc(N * sizeof(double)); // solution vector
    x = malloc(N * sizeof(double)); // variable vector

    // Coefficients for the tridiagonal matrix
    for (int i = 0; i < N; i++)
    {
        a[i] = -alpha;
        b[i] = 1 + 2 * alpha;
        c[i] = -alpha;
    }
    a[0] = 0;
    c[N - 1] = 0;

    // Check time
    double start, finish, elapsed;
    start = omp_get_wtime();

    // Shared variables
    int n;

    // Operator splitting method
    for (n = 0; n < M - 1; n++)
    {
        int n_exp, i;

        // Explicit method part: ODEs for u, v, w, s and other reaction functions
        for (n_exp = 0; n_exp < M_exp; n_exp++)
        {
            # pragma omp parallel for num_threads(num_threads) default(none) \
            private(i, I_app, tau_vminus, tau_wminus, tau_so, tau_s, tau_o, J_fi, J_so, J_si, J, v_inf, w_inf, du_dt, \
            dv_dt, ds_dt, dw_dt ) \
            shared(n, U, V, W, S, d, N, M, a, b, c, x, \
            u_o, u_u, theta_v, theta_w, theta_vminus, theta_o, tau_v1minus, tau_v2minus, \
            tau_vplus, tau_w1minus, tau_w2minus, k_wminus, u_wminus, tau_wplus, tau_fi, tau_o1, tau_o2, tau_so1, tau_so2, k_so, \
            u_so, tau_s1, tau_s2, k_s, u_s, tau_si, tau_winf, w_infstar, D, \
            M_exp, delta_t_exp, delta_t)

            for (i = 1; i < N - 1; i++)
            {
                // Stimulus
                if (n >= 3000 && n <= 5000 && i > 1 && i < 30)
                {
                    I_app = 0.5;
                }
                else
                {
                    I_app = 0;
                }

                tau_vminus = (1 - H(U[n * N + i], theta_vminus)) * tau_v1minus + H(U[n * N + i], theta_vminus) * tau_v2minus;
                tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1 + tanh(k_wminus * (U[n * N + i] - u_wminus))) / 2;
                tau_so = tau_so1 + (tau_so2 - tau_so1) * (1 + tanh(k_so * (U[n * N + i] - u_so))) / 2;
                tau_s = (1 - H(U[n * N + i], theta_w)) * tau_s1 + H(U[n * N + i], theta_w) * tau_s2;
                tau_o = (1 - H(U[n * N + i], theta_o)) * tau_o1 + H(U[n * N + i], theta_o) * tau_o2;

                J_fi = -V[i] * H(U[n * N + i], theta_v) * (U[n * N + i] - theta_v) * (u_u - U[n * N + i]) / tau_fi;
                J_so = (U[n * N + i] - u_o) * ((1 - H(U[n * N + i], theta_w)) / tau_o) + (H(U[n * N + i], theta_w) / tau_so);
                J_si = -H(U[n * N + i], theta_w) * W[i] * S[i] / tau_si;
                J = J_fi + J_so + J_si;

                v_inf = v_inf_function(U[n * N + i], theta_vminus);
                w_inf = (1 - H(U[n * N + i], theta_o)) * (1 - U[n * N + i] / tau_winf) + H(U[n * N + i], theta_o) * w_infstar;

                du_dt = -J + I_app;
                dv_dt = (1 - H(U[n * N + i], theta_v)) * (v_inf - V[i]) / tau_vminus - H(U[n * N + i], theta_v) * V[i] / tau_vplus;
                dw_dt = (1 - H(U[n * N + i], theta_w)) * (w_inf - W[i]) / tau_wminus - H(U[n * N + i], theta_w) * W[i] / tau_wplus;
                ds_dt = ((1 + tanh(k_s * (U[n * N + i] - u_s))) / 2 - S[i]) / tau_s;

                // Update array that will be used as input for the implicit method
                d[i] = U[n * N + i] + du_dt * delta_t_exp;

                // Update variables
                U[n * N + i] = U[n * N + i] + du_dt * delta_t_exp;
                V[i] = V[i] + dv_dt * delta_t_exp;
                W[i] = W[i] + dw_dt * delta_t_exp;
                S[i] = S[i] + ds_dt * delta_t_exp;
            }
        }

        d[0] = d[1];
        d[N - 1] = d[N - 2];

        // Implicit part: diffusion PDE (u)

        // Solve tridiagonal matrix (Linear system)
        x = solve_tridiagonal_in_place_reusable(d, N, a, b, c);

        // Update array U
        for (int k = 0; k < N; k++)
        {
            U[(n + 1) * N + k] = x[k];
        }

        // Boundary Condition
        /* U[(n + 1) * N] = U[(n + 1) * N + 1];
        U[(n + 1) * N + (N - 1)] = U[(n + 1) * N + (N - 2)]; */
    }

    // Check time
    finish = omp_get_wtime();
    elapsed = finish - start;
    printf("Elapsed time = %e seconds\n", elapsed);

    // Write to file simple
    FILE *fp = NULL;
    FILE *fp2 = NULL;

    fp = fopen("omp-minmodel4.txt", "w");
    fp2 = fopen("omp-minmodel10.txt", "w");
    for (int i = 0; i < N; i++)
    {
        fprintf(fp, "%lf\n", U[4000 * N + i]);
        fprintf(fp2, "%lf\n", U[10000 * N + i]);
    }
    printf("Files ready\n");

    FILE *fp_all = NULL;
    fp_all = fopen("omp-minmodel-all.txt", "w");
    int count = 0;
    bool tag = false;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i % 100 == 0)
            {
                fprintf(fp_all, "%lf\n", U[i * N + j]);
                tag = true;
            }
        }
        if (tag)
        {
            count++;
            tag = false;
        }
    }
    printf("File complete ready with c = %d\n", count);

    free(U);
    free(V);
    free(W);
    free(S);

    return 0;
}
