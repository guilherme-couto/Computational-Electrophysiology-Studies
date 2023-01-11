#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

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

int main(int argc, char *argv[])
{
    // Discretization
    int T = 600;
    double L = 300;
    double delta_x = 1;
    double delta_t = 0.001;

    int M = T / delta_t; // number of points in time
    int N = L / delta_x; // number of points in space

    printf("Number of points in time = %d\n", M);
    printf("Number of points in space = %d\n", N);

    // Parameters
    double *U, *V, *W, *S;
    U = malloc(M * N * sizeof(double));
    V = malloc(M * N * sizeof(double));
    W = malloc(M * N * sizeof(double));
    S = malloc(M * N * sizeof(double));

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

    // Explicit Euler method
    for (int n = 0; n < M - 1; n++)
    {
        for (int i = 1; i < N - 1; i++)
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

            J_fi = -V[n * N + i] * H(U[n * N + i], theta_v) * (U[n * N + i] - theta_v) * (u_u - U[n * N + i]) / tau_fi;
            J_so = (U[n * N + i] - u_o) * ((1 - H(U[n * N + i], theta_w)) / tau_o) + (H(U[n * N + i], theta_w) / tau_so);
            J_si = -H(U[n * N + i], theta_w) * W[n * N + i] * S[n * N + i] / tau_si;
            J = J_fi + J_so + J_si;

            v_inf = v_inf_function(U[n * N + i], theta_vminus);
            w_inf = (1 - H(U[n * N + i], theta_o)) * (1 - U[n * N + i] / tau_winf) + H(U[n * N + i], theta_o) * w_infstar;

            du_dt = -J + I_app;
            dv_dt = (1 - H(U[n * N + i], theta_v)) * (v_inf - V[n * N + i]) / tau_vminus - H(U[n * N + i], theta_v) * V[n * N + i] / tau_vplus;
            dw_dt = (1 - H(U[n * N + i], theta_w)) * (w_inf - W[n * N + i]) / tau_wminus - H(U[n * N + i], theta_w) * W[n * N + i] / tau_wplus;
            ds_dt = ((1 + tanh(k_s * (U[n * N + i] - u_s))) / 2 - S[n * N + i]) / tau_s;

            // Update variables
            U[(n + 1) * N + i] = U[n * N + i] + (D * (U[n * N + (i + 1)] - 2 * U[n * N + i] + U[n * N + (i - 1)]) / delta_x / delta_x + du_dt) * delta_t;
            V[(n + 1) * N + i] = V[n * N + i] + dv_dt * delta_t;
            W[(n + 1) * N + i] = W[n * N + i] + dw_dt * delta_t;
            S[(n + 1) * N + i] = S[n * N + i] + ds_dt * delta_t;
        }

        // Boundary Condition
        U[(n + 1) * N] = U[(n + 1) * N + 1];
        U[(n + 1) * N + (N - 1)] = U[(n + 1) * N + (N - 2)];
    }

    // Write to file simple
    FILE *fp = NULL;
    FILE *fp2 = NULL;

    fp = fopen("ceq-minmodel4.txt", "w");
    fp2 = fopen("ceq-minmodel10.txt", "w");
    for (int i = 0; i < N; i++)
    {
        fprintf(fp, "%lf\n", U[4000 * N + i]);
        fprintf(fp2, "%lf\n", U[10000 * N + i]);
    }
    printf("Files ready\n");

    FILE *fp_all = NULL;
    fp_all = fopen("ceq-minmodel-all.txt", "w");
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
