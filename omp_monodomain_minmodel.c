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

void thomas_algorithm(double *d, unsigned long N, double alpha)
{
    // Auxiliary arrays
    double *c_ = (double *)malloc((N-3) * sizeof(double));
    double *d_ = (double *)malloc((N-2) * sizeof(double));

    // Coefficients
    double a = -alpha;          // subdiagonal
    double b = 1 + alpha;       // diagonal (1st and last row)
    double c = -alpha;          // superdiagonal

    // 1st: update auxiliary arrays
    c_[0] = c / b;
    d_[0] = d[1] / b;

    b = 1 + 2*alpha;

    for (int i = 1; i < N-3; i++)
    {
        c_[i] = c / (b - a * c_[i-1]);
        d_[i] = (d[i+1] - a * d_[i-1]) / (b - a * c_[i-1]);
    }

    b = 1 + alpha;
    d_[N-3] = (d[N-2] - a * d_[N-4]) / (b - a * c_[N-4]);

    // 2nd: update solution
    d[N-2] = d_[N-3];

    for (int i = N-3; i >= 1; i--)
    {
        d[i] = d_[i-1] - c_[i-1] * d[i+1];
    }

    // Free memory
    free(c_);
    free(d_);
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
    int L = 150;
    double delta_x = 1;
    double delta_y = 1;
    double delta_t = 0.02;
    double delta_t_ode = 0.01;

    int M = T / delta_t;                // number of points in time
    int M_ode = delta_t / delta_t_ode;  // number of points in time for ODE
    int N_x = L / delta_x;              // number of points in x
    int N_y = L / delta_y;              // number of points in y

    printf("Number of points in time = %d\n", M);
    printf("ODE iterations for each point in time = %d\n", M_ode);
    printf("Number of points in x = %d\n", N_x);
    printf("Number of points in y = %d\n\n", N_y);

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

    // Ask for input
    char method;
    if (delta_t == delta_t_ode)
    {
        printf("Which method (e - explicit | a - ADI): ");
        scanf("%c", &method);
    }
    else
        method = 'a';

    // Shared variables
    int n, n_ode;
    int t_app = 2 / delta_t;

    // Open the file to write
    /* FILE *fp_all = NULL;
    fp_all = fopen("omp-mono-all.txt", "w");
    int count = 0; */

    FILE *fp = NULL, *fp2 = NULL;
    fp = fopen("mono_test_exp.txt", "w");
    fp2 = fopen("mono_test_adi.txt", "w");
    //fp = fopen("minmodel_times.txt", "a");

    // Start timer
    double start, finish, elapsed;
    start = omp_get_wtime();

    if (method == 'e' || method == 'E')
    {
        // Explicit (parallelized)
        for (n = 0; n < M - 1; n++)
        {
            int i, j;

            # pragma omp parallel for collapse(2) num_threads(num_threads) default(none) \
            private(i, j, I_app, tau_vminus, tau_wminus, tau_so, tau_s, tau_o, J_fi, J_so, J_si, J, v_inf, w_inf, du_dt, \
            dv_dt, ds_dt, dw_dt ) \
            shared(t_app, n, n_ode, U, V, W, S, N_x, N_y, M, M_ode, \
            u_o, u_u, theta_v, theta_w, theta_vminus, theta_o, tau_v1minus, tau_v2minus, \
            tau_vplus, tau_w1minus, tau_w2minus, k_wminus, u_wminus, tau_wplus, tau_fi, tau_o1, tau_o2, tau_so1, tau_so2, k_so, \
            u_so, tau_s1, tau_s2, k_s, u_s, tau_si, tau_winf, w_infstar, D, \
            delta_t, delta_x, delta_y)

            for (i = 1; i < N_x - 1; i++)
            {
                for (j = 1; j < N_y - 1; j++)
                {
                    // Calculate reaction functions (ODE)

                    // Stimulus
                    if (n >= 0 && n <= t_app && j > 0 && j < 10)
                        I_app = 1;
                    else
                        I_app = 0;

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

                    // Update variables
                    U[i][j] = U[i][j] + du_dt * delta_t;
                    V[i][j] = V[i][j] + dv_dt * delta_t;
                    W[i][j] = W[i][j] + dw_dt * delta_t;
                    S[i][j] = S[i][j] + ds_dt * delta_t;
                    
                    // Diffusion (PDE)
                    U[i][j] = U[i][j] + (delta_t * D) * ((U[i - 1][j] - 2 * U[i][j] + U[i + 1][j]) / delta_x / delta_x);
                    U[i][j] = U[i][j] + (delta_t * D) * ((U[i][j - 1] - 2 * U[i][j] + U[i][j + 1]) / delta_y / delta_y);
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

            // Write to file
            /* if (n % 10 == 0)
            {
                for (int i = 0; i < N_x; i++)
                {
                    for (int j = 0; j < N_y; j++)
                    {
                        fprintf(fp_all, "%lf\n", U[i][j]);
                    }
                }
                count++;
            } */

            // Write to file
            // Error analysis
            if (n * delta_t == 42)
            {
                for (int i = 0; i < N_x; i++)
                {
                    for (int j = 0; j < N_y; j++)
                    {
                        fprintf(fp, "%lf\n", U[i][j]);
                    }
                    
                }
            }
        }
    }
    else if (method == 'a' || method == 'A')
    {
        // ADI Method
        for (n = 0; n < M - 1; n++)
        {
            int i, j;

            // STEP 1: implicit in x, explicit in y on interval [t_n, t_n+(1/2)]

            # pragma omp parallel for collapse(2) num_threads(num_threads) default(none) \
            private(i, j, I_app, tau_vminus, tau_wminus, tau_so, tau_s, tau_o, J_fi, J_so, J_si, J, v_inf, w_inf, du_dt, \
            dv_dt, ds_dt, dw_dt, n_ode) \
            shared(n, U, V, W, S, N_x, N_y, M_ode, \
            u_o, u_u, theta_v, theta_w, theta_vminus, theta_o, tau_v1minus, tau_v2minus, \
            tau_vplus, tau_w1minus, tau_w2minus, k_wminus, u_wminus, tau_wplus, tau_fi, tau_o1, tau_o2, tau_so1, tau_so2, k_so, \
            u_so, tau_s1, tau_s2, k_s, u_s, tau_si, tau_winf, w_infstar, \
            delta_t, alpha, delta_t_ode, num_threads, t_app)
            for (i = 1; i < N_x - 1; i++)
            {
                for (j = 1; j < N_y - 1; j++)
                {
                    // ODEs
                    for (n_ode = 0; n_ode < M_ode; n_ode++)
                    {
                        // Stimulus
                        if (n >= 0 && n <= t_app && j > 0 && j < 10)
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
                        U[i][j] = U[i][j] + delta_t_ode * du_dt;

                        // Update gating variables (v, w, s)
                        V[i][j] = V[i][j] + delta_t_ode * dv_dt;
                        W[i][j] = W[i][j] + delta_t_ode * dw_dt;
                        S[i][j] = S[i][j] + delta_t_ode * ds_dt;
                    }
                    U[i][j] = U[i][j] + (alpha * (U[i][j - 1] - 2*U[i][j] + U[i][j + 1]));
                }
            }

            // Solve tridiagonal matrix (Linear system) for x-axis
            # pragma omp parallel for num_threads(num_threads) default(none) \
            private(i) \
            shared(N_x, N_y, U, alpha)
            for (i = 1; i < N_x-1; i++)
                thomas_algorithm(U[i], N_y, alpha);
            
            // Boundary Conditions x-axis
            # pragma omp parallel for num_threads(num_threads) default(none) \
            private(j) \
            shared(U, N_x, N_y)
            for (j = 0; j < N_y; j++)
            {
                U[0][j] = U[1][j];
                U[N_x - 1][j] = U[N_x - 2][j];
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

            // STEP 2: explicit in x, implicit in y on interval [t_n+(1/2), t_n+1]

            # pragma omp parallel for collapse(2) num_threads(num_threads) default(none) \
            private(i, j) \
            shared(U, N_x, N_y, alpha)
            for (i = 1; i < N_x - 1; i++)
            {
                for (j = 1; j < N_y - 1; j++)
                {
                    // Potential (u)
                    U[i][j] = U[i][j] + (alpha * (U[i-1][j] - 2*U[i][j] + U[i+1][j]));
                }
            }

            // Solve tridiagonal matrix (Linear system) for y-axis
            # pragma omp parallel for num_threads(num_threads) default(none) \
            private(i) \
            shared(U, N_x, N_y, alpha)
            for (i = 1; i < N_x-1; i++)
                thomas_algorithm(U[i], N_y, alpha);
            
            // Boundary Conditions x-axis
            # pragma omp parallel for num_threads(num_threads) default(none) \
            private(j) \
            shared(U, N_x, N_y)
            for (j = 0; j < N_y; j++)
            {
                U[0][j] = U[1][j];
                U[N_x - 1][j] = U[N_x - 2][j];
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

            // Error analysis
            // Write to file
            if (n*delta_t == 42)
            {
                for (int i = 0; i < N_x; i++)
                {
                    for (int j = 0; j < N_y; j++)
                    {
                        fprintf(fp2, "%lf\n", U[i][j]);
                    }
                }
            }
            
            // Write to file
            /* if (n % 10 == 0)
            {
                for (int i = 0; i < N_x; i++)
                {
                    for (int j = 0; j < N_y; j++)
                    {
                        fprintf(fp_all, "%lf\n", U[i][j]);
                    }
                }
                count++;
            } */
        }
    }

    // Check time
    finish = omp_get_wtime();
    elapsed = finish - start;

    //fprintf(fp, "%e\n", elapsed);

    printf("\nElapsed time = %e seconds\n", elapsed);
    //printf("File complete ready with time dimension c = %d\n", count);

    fclose(fp);
    fclose(fp2);

    free(U);
    free(V);
    free(W);
    free(S);

    return 0;
}
