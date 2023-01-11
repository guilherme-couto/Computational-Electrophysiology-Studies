#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

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
    double delta_t = 0.001;

    int M = T / delta_t; // number of points in time

    printf("Number of points in time = %d\n", M);

    // Parameters
    double *U, *V, *W, *S;
    U = malloc(M * sizeof(double));
    V = malloc(M * sizeof(double));
    W = malloc(M * sizeof(double));
    S = malloc(M * sizeof(double));

    double J_fi, J_so, J_si, J, I_app;
    double tau_vminus, tau_wminus, tau_so, tau_s, tau_o;
    double v_inf, w_inf;
    double du_dt, dv_dt, ds_dt, dw_dt;

    // Initial conditions
    U[0] = 0;
    V[0] = 1;
    W[0] = 1;
    S[0] = 0;

    // Explicit Euler method
    for (int n = 0; n < M - 1; n++)
    {
        // Stimulus
        if (n >= 3000 && n <= 5000)
        {
            I_app = 0.5;
        }
        else
        {
            I_app = 0;
        }

        tau_vminus = (1 - H(U[n], theta_vminus)) * tau_v1minus + H(U[n], theta_vminus) * tau_v2minus;
        tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1 + tanh(k_wminus * (U[n] - u_wminus))) / 2;
        tau_so = tau_so1 + (tau_so2 - tau_so1) * (1 + tanh(k_so * (U[n] - u_so))) / 2;
        tau_s = (1 - H(U[n], theta_w)) * tau_s1 + H(U[n], theta_w) * tau_s2;
        tau_o = (1 - H(U[n], theta_o)) * tau_o1 + H(U[n], theta_o) * tau_o2;

        J_fi = -V[n] * H(U[n], theta_v) * (U[n] - theta_v) * (u_u - U[n]) / tau_fi;
        J_so = (U[n] - u_o) * ((1 - H(U[n], theta_w)) / tau_o) + (H(U[n], theta_w) / tau_so);
        J_si = -H(U[n], theta_w) * W[n] * S[n] / tau_si;
        J = J_fi + J_so + J_si;

        v_inf = v_inf_function(U[n], theta_vminus);
        w_inf = (1 - H(U[n], theta_o)) * (1 - U[n] / tau_winf) + H(U[n], theta_o) * w_infstar;

        du_dt = -J + I_app;
        dv_dt = (1 - H(U[n], theta_v)) * (v_inf - V[n]) / tau_vminus - H(U[n], theta_v) * V[n] / tau_vplus;
        dw_dt = (1 - H(U[n], theta_w)) * (w_inf - W[n]) / tau_wminus - H(U[n], theta_w) * W[n] / tau_wplus;
        ds_dt = ((1 + tanh(k_s * (U[n] - u_s))) / 2 - S[n]) / tau_s;

        // Update variables
        U[n + 1] = U[n] + du_dt * delta_t;
        V[n + 1] = V[n] + dv_dt * delta_t;
        W[n + 1] = W[n] + dw_dt * delta_t;
        S[n + 1] = S[n] + ds_dt * delta_t;
    }

    // Boundary Condition
    U[0] = U[1];
    U[M - 1] = U[M - 2];

    // Write to file
    FILE *fp = NULL;
    fp = fopen("minmodel.txt", "w");
    double t = 0;
    for (int i = 0; i < M; i++)
    {
        fprintf(fp, "%lf\n", U[i]);
        t += delta_t;
    }
    printf("File ready\n");

    free(U);
    free(V);
    free(W);
    free(S);

    return 0;
}
