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
double beta = 1;  // surface to volume ratio

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

void thomas(double *d, unsigned long N, double alpha)
{
    // Constants
    double a = -alpha;
    double b = (1 + 2 * alpha);
    double c = - alpha;

    c = c / b;
    d[0] = d[0] / b;

    /* loop from 1 to N - 1 inclusive */
    for (int ix = 1; ix < N; ix++)
    {
        if (ix == N - 1)
           c = 0;
        else
            c = -alpha;

        double m = 1.0 / (b - a * c);
        
        d[ix] = (d[ix] - a * d[ix - 1]) * m;
    }

    /* loop from N - 2 to 0 inclusive, safely testing loop end condition */
    for (int ix = N - 2; ix >= 0; ix--)
    {
        if (ix == 0)
            a = 0;
        else
            a = -alpha;

        if (ix == N - 1)
           c = 0;
        else
            c = -alpha;
        
        double m = 1.0 / (b - a * c);
        c = c * m;

        d[ix] = d[ix] - (c * d[ix + 1]);
    }
    return;
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
    int T = 100;
    int L = 100;
    double delta_x = 1;
    double delta_y = 1;
    double delta_t = 0.001;

    int M = T / delta_t;               // number of points in time
    int N_x = L / delta_x;             // number of points in x
    int N_y = L / delta_y;             // number of points in y

    printf("Number of points in time = %d\n", M);
    printf("Number of points in x = %d\n", N_x);
    printf("Number of points in y = %d\n", N_y);

    // Parameters
    double **U = (double **)malloc(N_x * sizeof(double *));
    double **V = (double **)malloc(N_x * sizeof(double *));
    double **W = (double **)malloc(N_x * sizeof(double *));
    double **S = (double **)malloc(N_x * sizeof(double *));

    for (int i = 0; i < N_x; i++)
    {
        U[i] = (double *)malloc(N_y * sizeof(double));
        V[i] = (double *)malloc(N_y * sizeof(double));
        W[i] = (double *)malloc(N_y * sizeof(double));
        S[i] = (double *)malloc(N_y * sizeof(double));
    }

    double J_fi, J_so, J_si, J = 0, I_app;
    double tau_vminus, tau_wminus, tau_so, tau_s, tau_o;
    double v_inf, w_inf;
    double du_dt = 0, dv_dt = 0, ds_dt = 0, dw_dt = 0;

    // Initial conditions
    for (int i = 0; i < N_x; i++)
    {
        for (int j = 0; j < N_y; j++)
        {
            U[i][j] = 0; // i: x-axis, j: y-axis
            V[i][j] = 1;
            W[i][j] = 1;
            S[i][j] = 0;
        }
    }

    // Subdiagonal, diagonal and superdiagonal of the tridiagonal matrix (Implicit Method)
    double alpha = (D * delta_t) / (2 * delta_x * delta_x);

    // Check time
    double start, finish, elapsed;
    start = omp_get_wtime();

    // Shared variables
    int n;

    // Write to file simple
    FILE *fp_all = NULL;
    fp_all = fopen("omp-mono-all.txt", "w");
    int count = 0;

    // Explicit with Splitting (parallelized)
    /* for (n = 0; n < M - 1; n++)
    {
        int i, k;

        // Calculate reaction functions
        # pragma omp parallel for collapse(2) num_threads(num_threads) default(none) \
        private(i, k, I_app, tau_vminus, tau_wminus, tau_so, tau_s, tau_o, J_fi, J_so, J_si, J, v_inf, w_inf, du_dt, \
        dv_dt, ds_dt, dw_dt ) \
        shared(n, U, V, W, S, N_x, N_y, M, \
        u_o, u_u, theta_v, theta_w, theta_vminus, theta_o, tau_v1minus, tau_v2minus, \
        tau_vplus, tau_w1minus, tau_w2minus, k_wminus, u_wminus, tau_wplus, tau_fi, tau_o1, tau_o2, tau_so1, tau_so2, k_so, \
        u_so, tau_s1, tau_s2, k_s, u_s, tau_si, tau_winf, w_infstar, D, beta, \
        delta_t, delta_x, delta_y)

        for (i = 1; i < N_x - 1; i++)
        {
            for (k = 1; k < N_y - 1; k++)
            {
                // Stimulus
                if (n >= 0 && n <= 2000 && k > 0 && k < 10)
                    I_app = 1;
                else
                    I_app = 0;

                tau_vminus = (1 - H(U[i][k], theta_vminus)) * tau_v1minus + H(U[i][k], theta_vminus) * tau_v2minus;
                tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1 + tanh(k_wminus * (U[i][k] - u_wminus))) / 2;
                tau_so = tau_so1 + (tau_so2 - tau_so1) * (1 + tanh(k_so * (U[i][k] - u_so))) / 2;
                tau_s = (1 - H(U[i][k], theta_w)) * tau_s1 + H(U[i][k], theta_w) * tau_s2;
                tau_o = (1 - H(U[i][k], theta_o)) * tau_o1 + H(U[i][k], theta_o) * tau_o2;

                J_fi = -V[i][k] * H(U[i][k], theta_v) * (U[i][k] - theta_v) * (u_u - U[i][k]) / tau_fi;
                J_so = (U[i][k] - u_o) * ((1 - H(U[i][k], theta_w)) / tau_o) + (H(U[i][k], theta_w) / tau_so);
                J_si = -H(U[i][k], theta_w) * W[i][k] * S[i][k] / tau_si;
                J = J_fi + J_so + J_si;

                v_inf = v_inf_function(U[i][k], theta_vminus);
                w_inf = (1 - H(U[i][k], theta_o)) * (1 - U[i][k] / tau_winf) + H(U[i][k], theta_o) * w_infstar;

                du_dt = -J + I_app;
                dv_dt = (1 - H(U[i][k], theta_v)) * (v_inf - V[i][k]) / tau_vminus - H(U[i][k], theta_v) * V[i][k] / tau_vplus;
                dw_dt = (1 - H(U[i][k], theta_w)) * (w_inf - W[i][k]) / tau_wminus - H(U[i][k], theta_w) * W[i][k] / tau_wplus;
                ds_dt = ((1 + tanh(k_s * (U[i][k] - u_s))) / 2 - S[i][k]) / tau_s;

                // Update variables
                U[i][k] = U[i][k] + (du_dt * delta_t);
                V[i][k] = V[i][k] + dv_dt * delta_t;
                W[i][k] = W[i][k] + dw_dt * delta_t;
                S[i][k] = S[i][k] + ds_dt * delta_t;

                // Diffusion
                // x-axis
                U[i][k] = U[i][k] + (delta_t * (D / beta) * ((U[i - 1][k] - 2 * U[i][k] + U[i + 1][k]) / delta_x / delta_x));

                // y-axis
                U[i][k] = U[i][k] + (delta_t * (D / beta) * ((U[i][k - 1] - 2 * U[i][k] + U[i][k + 1]) / delta_y / delta_y));
            }
        }

        // Boundary Conditions y-axis
        # pragma omp parallel for num_threads(num_threads) default(none) \
        private(i) \
        shared(U, N_x, N_y)

        for (i = 0; i < N_x; i++)
        {
            U[i][0] = U[i][1];
            U[i][N_y - 1] = U[i][N_y - 2];
        }

        // Boundary Conditions x-axis
        # pragma omp parallel for num_threads(num_threads) default(none) \
        private(k) \
        shared(U, N_x, N_y)

        for (k = 0; k < N_y; k++)
        {
            U[0][k] = U[1][k];
            U[N_x - 1][k] = U[N_x - 2][k];
        }

        // Write to file
        if (n % 10 == 0)
        {
            for (int i = 0; i < N_x; i++)
            {
                for (int k = 0; k < N_y; k++)
                {
                    fprintf(fp_all, "%lf\n", U[i][k]);
                }
            }
            count++;
        }
    } */

    // ADI Method
    for (n = 0; n < M - 1; n++)
    {
        int i, j;

        // STEP 1: implicit in x, explicit in y on interval [t_n, t_n+(1/2)]

        # pragma omp parallel for collapse(2) num_threads(num_threads) default(none) \
        private(i, j, I_app, tau_vminus, tau_wminus, tau_so, tau_s, tau_o, J_fi, J_so, J_si, J, v_inf, w_inf, du_dt, \
        dv_dt, ds_dt, dw_dt ) \
        shared(n, U, V, W, S, N_x, N_y, \
        u_o, u_u, theta_v, theta_w, theta_vminus, theta_o, tau_v1minus, tau_v2minus, \
        tau_vplus, tau_w1minus, tau_w2minus, k_wminus, u_wminus, tau_wplus, tau_fi, tau_o1, tau_o2, tau_so1, tau_so2, k_so, \
        u_so, tau_s1, tau_s2, k_s, u_s, tau_si, tau_winf, w_infstar, beta, \
        delta_t, alpha)
        for (i = 1; i < N_x - 1; i++)
        {
            for (j = 1; j < N_y - 1; j++)
            {
                // Stimulus
                if (n >= 0 && n <= 2000 && j > 0 && j < 10)
                    I_app = 1;
                else
                    I_app = 0;

                // Calculate J (I_ion) functions
                tau_vminus = (1 - H(U[i][j], theta_vminus)) * tau_v1minus + H(U[i][j], theta_vminus) * tau_v2minus;
                tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1 + tanh(k_wminus * (U[i][j] - u_wminus))) / 2;
                tau_so = tau_so1 + (tau_so2 - tau_so1) * (1 + tanh(k_so * (U[i][j] - u_so))) / 2;
                tau_s = (1 - H(U[i][j], theta_w)) * tau_s1 + H(U[i][j], theta_w) * tau_s2;
                tau_o = (1 - H(U[i][j], theta_o)) * tau_o1 + H(U[i][j], theta_o) * tau_o2;

                J_fi = -V[i][j] * H(U[i][j], theta_v) * (U[i][j] - theta_v) * (u_u - U[i][j]) / tau_fi;
                J_so = (U[i][j] - u_o) * ((1 - H(U[i][j], theta_w)) / tau_o) + (H(U[i][j], theta_w) / tau_so);
                J_si = -H(U[i][j], theta_w) * W[i][j] * S[i][j] / tau_si;
                J = J_fi + J_so + J_si;

                v_inf = v_inf_function(U[i][j], theta_vminus);
                w_inf = (1 - H(U[i][j], theta_o)) * (1 - U[i][j] / tau_winf) + H(U[i][j], theta_o) * w_infstar;

                du_dt = -J + I_app;
                dv_dt = (1 - H(U[i][j], theta_v)) * (v_inf - V[i][j]) / tau_vminus - H(U[i][j], theta_v) * V[i][j] / tau_vplus;
                dw_dt = (1 - H(U[i][j], theta_w)) * (w_inf - W[i][j]) / tau_wminus - H(U[i][j], theta_w) * W[i][j] / tau_wplus;
                ds_dt = ((1 + tanh(k_s * (U[i][j] - u_s))) / 2 - S[i][j]) / tau_s;

                // Potential (u)
                U[i][j] = U[i][j] + (du_dt * delta_t / 2) + (alpha * (U[i][j - 1] - 2*U[i][j] + U[i][j + 1]));

                // Update gating variables (v, w, s)
                V[i][j] = V[i][j] + delta_t * dv_dt;
                W[i][j] = W[i][j] + delta_t * dw_dt;
                S[i][j] = S[i][j] + delta_t * ds_dt;
            }
        }

        // Boundary Conditions y-axis
        # pragma omp parallel for num_threads(num_threads) default(none) \
        private(i) \
        shared(U, N_x, N_y)
        for (i = 0; i < N_x; i++)
        {
            U[i][0] = U[i][1];
            U[i][N_y - 1] = U[i][N_y - 2];
        }

        // Boundary Conditions x-axis
        # pragma omp parallel for num_threads(num_threads) default(none) \
        private(j) \
        shared(U, N_x, N_y)
        for (j = 0; j < N_y; j++)
        {
            U[0][j] = U[1][j];
            U[N_x - 1][j] = U[N_x - 2][j];
        }

        // Solve tridiagonal matrix (Linear system) for x-axis
        # pragma omp parallel for num_threads(num_threads) default(none) \
        private(i) \
        shared(N_x, N_y, U, alpha)
        for (i = 0; i < N_x; i++)
            thomas(U[i], N_y, alpha);


        // STEP 2: explicit in x, implicit in y on interval [t_n+(1/2), t_n+1]

        # pragma omp parallel for collapse(2) num_threads(num_threads) default(none) \
        private(i, j) \
        shared(U, N_x, N_y, delta_t, alpha)
        for (i = 1; i < N_x - 1; i++)
        {
            for (j = 1; j < N_y - 1; j++)
            {
                // Potential (u)
                U[i][j] = U[i][j] + (alpha * (U[i-1][j] - 2*U[i][j] + U[i+1][j]));
            }
        }

        // Boundary Conditions y-axis
        # pragma omp parallel for num_threads(num_threads) default(none) \
        private(i) \
        shared(U, N_x, N_y)
        for (i = 0; i < N_x; i++)
        {
            U[i][0] = U[i][1];
            U[i][N_y - 1] = U[i][N_y - 2];
        }

        // Boundary Conditions x-axis
        # pragma omp parallel for num_threads(num_threads) default(none) \
        private(j) \
        shared(U, N_x, N_y)
        for (j = 0; j < N_y; j++)
        {
            U[0][j] = U[1][j];
            U[N_x - 1][j] = U[N_x - 2][j];
        }

        // Solve tridiagonal matrix (Linear system) for y-axis
        # pragma omp parallel for num_threads(num_threads) default(none) \
        private(i) \
        shared(U, N_x, N_y, alpha)
        for (i = 0; i < N_x; i++)
            thomas(U[i], N_y, alpha);


        // Write to file
        if (n % 10 == 0)
        {
            for (int i = 0; i < N_x; i++)
            {
                for (int j = 0; j < N_y; j++)
                {
                    fprintf(fp_all, "%lf\n", U[i][j]);
                }
            }
            count++;
        }
    }

    // Check time
    finish = omp_get_wtime();
    elapsed = finish - start;
    printf("Elapsed time = %e seconds\n", elapsed);

    printf("File complete ready with time dimension c = %d\n", count);

    free(U);
    free(V);
    free(W);
    free(S);

    return 0;
}
