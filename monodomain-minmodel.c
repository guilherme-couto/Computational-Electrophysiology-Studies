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

double D = 1.171;   // +- 0.0221 cm^2/s  human ventricular diffusion coefficient

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
void solve_tridiagonal_in_place_reusable(double * x, int X, double * a, double * b, double * c) {
    /* Allocate scratch space. */
    double * cprime = malloc(sizeof(double) * X);
    
    if (!cprime) {
        printf("Error: malloc failed.\n");
    }
    cprime[0] = c[0] / b[0];
    x[0] = x[0] / b[0];
    
    /* loop from 1 to X - 1 inclusive */
    for (int ix = 1; ix < X; ix++) {
        double m = 1.0 / (b[ix] - a[ix] * cprime[ix - 1]);
        cprime[ix] = c[ix] * m;
        x[ix] = (x[ix] - a[ix] * x[ix - 1]) * m;
    }
    
    /* loop from X - 2 to 0 inclusive, safely testing loop end condition */
    for (int ix = X - 2; ix >= 0; ix--)
        x[ix] = x[ix] - cprime[ix] * x[ix + 1];
    
    /* free scratch space */
    free(cprime);
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
    int T = 1;
    double L = 1200;
    double delta_x = 1;
    double delta_t = 0.0001;
    double delta_t_exp = 0.0001;

    int M = T / delta_t;                // number of points in time
    int N = L / delta_x;                // number of points in space
    int M_exp = delta_t / delta_t_exp;  // number of points in space for explicit method

    printf("Number of points in time = %d\n", M);
    printf("Number of points in space = %d\n", N);
    printf("Number of points in time for explicit method = %d\n", M_exp);

    // Parameters
    double *U, *V, *W, *S;
    U = malloc(M * N * sizeof(double));
    V = malloc(M * N * sizeof(double));
    W = malloc(M * N * sizeof(double));
    S = malloc(M * N * sizeof(double));

    //double u, v, w, s;
    double J_fi, J_so, J_si, J;
    double tau_vminus, tau_wminus, tau_so, tau_s, tau_o;
    double v_inf, w_inf;
    double du_dt, dv_dt, ds_dt, dw_dt;

    // Initial conditions
    // u = 0;
    // v = 1;
    // w = 1;
    // s = 0;

    for (int i = 0; i < N; i++)
    {
        U[i] = 0;
        V[i] = 1;
        W[i] = 1;
        S[i] = 0;
    }

    // Subdiagonal, diagonal and superdiagonal of the tridiagonal matrix (Implicit Method)
    double alpha = (D*delta_t) / (delta_x * delta_x);
    double *a, *b, *c, *x;
    a = malloc(N * sizeof(double));     // subdiagonal
    b = malloc(N * sizeof(double));     // diagonal
    c = malloc(N * sizeof(double));     // superdiagonal
    x = malloc(N * sizeof(double));     // both variable and solution vector

    // Coefficients for the tridiagonal matrix
    for(int i = 0; i < N; i++)
    {
        a[i] = alpha;
        b[i] = 1 - 2*alpha;
        c[i] = alpha;
    }
    a[0] = 0;
    c[N-1] = 0;

    // Operator splitting method
    for (int n = 0; n < M - 1; n++)
    { 
        int n_exp, i;

        int step_n = n*N;
        int step_np1 = (n+1)*N;

        if(n == 1000)
        {
            U[n*N+1] = 0.325;
        } 

        // Explicit method part: ODEs for u, v, w, s and other reaction functions
        for (n_exp = 0; n_exp < M_exp; n_exp++)
        {
            for(i = 1; i < N-1; i++)
            {
                tau_vminus = (1-H(U[step_n+i], theta_vminus))*tau_v1minus + H(U[step_n+i], theta_vminus)*tau_v2minus;
                tau_wminus = tau_w1minus + (tau_w2minus-tau_w1minus)*(1+tanh(k_wminus*(U[step_n+i]-u_wminus)))/2;
                tau_so = tau_so1 + (tau_so2-tau_so1)*(1+tanh(k_so*(U[step_n+i]-u_so)))/2;
                tau_s = (1-H(U[step_n+i], theta_w))*tau_s1 + H(U[step_n+i], theta_w)*tau_s2;
                tau_o = (1-H(U[step_n+i], theta_o))*tau_o1 + H(U[step_n+i], theta_o)*tau_o2;

                J_fi = -V[step_n+i] * H(U[step_n+i], theta_v) * (U[step_n+i]-theta_v) * (u_u-U[step_n+i]) / tau_fi;
                J_so = (U[step_n+i]-u_o) * ((1-H(U[step_n+i], theta_w)) / tau_o) + (H(U[step_n+i], theta_w) / tau_so);
                J_si = -H(U[step_n+i], theta_w) * W[step_n+i] * S[step_n+i] / tau_si;
                J = J_fi + J_so + J_si;
                
                v_inf = v_inf_function(U[step_n+i], theta_vminus);
                w_inf = (1-H(U[step_n+i], theta_o))*(1-U[step_n+i]/tau_winf) + H(U[step_n+i], theta_o)*w_infstar;
                
                du_dt = - J;
                dv_dt = (1-H(U[step_n+i], theta_v))*(v_inf-V[step_n+i])/tau_vminus - H(U[step_n+i], theta_v)*V[step_n+i]/tau_vplus;
                dw_dt = (1-H(U[step_n+i], theta_w))*(w_inf-W[step_n+i])/tau_wminus - H(U[step_n+i], theta_w)*W[step_n+i]/tau_wplus;
                ds_dt = ((1+tanh(k_s*(U[step_n+i]-u_s)))/2 - S[step_n+i])/tau_s;
                
                // Update array that will be used as input for the implicit method
                x[i] = U[step_n+i];

                // Update variables
                U[step_np1+i] = U[step_n+i] + du_dt * delta_t_exp;
                V[step_np1+i] = V[step_n+i] + dv_dt * delta_t_exp;
                W[step_np1+i] = W[step_n+i] + dw_dt * delta_t_exp;
                S[step_np1+i] = S[step_n+i] + ds_dt * delta_t_exp;
            }       
        }
        x[0] = x[1];
        x[N-1] = x[N-2];

        // Implicit part: diffusion PDE (u)

        // Solve tridiagonal matrix (Linear system)
        solve_tridiagonal_in_place_reusable(x, N, a, b, c);
        
        // Update array U
        
        for(int i = 0; i < N; i++)
        {
            if(x[i] > 0.0){
                printf("%lf   ", x[i]);
            }
            U[step_np1 + i] = x[i];
        }

        // Boundary Condition
        U[step_np1 + 0] = U[step_np1 + 1];
        U[step_np1 + (N - 1)] = U[step_np1 + (N - 2)];

    } 


    FILE *fp = NULL;
    fp = fopen("ceq-mv.dat", "w");
    double t = 0;
    for (int i = 0; i < N; i++)
    {
        fprintf(fp, "%lf\t %lf\n", t, U[1000*N+i]);    
        t += delta_t;
    }
    printf("File ready\n");

    free(U);
    free(V);
    free(W);
    free(S);

    free(a);
    free(b);
    free(c);
    free(x);

    return 0;
}
