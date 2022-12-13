#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

// Parameters to ajust for TNNP model
double u_o = 0;
double u_u = 1.58;
double theta_v = 0.3;
double theta_w = 0.015;
double theta_vminus = 0.015;
double theta_o = 0.006;
double tau_v1minus = 60;
double tau_v2minus = 1150;
double tau_vplus = 1.4506;
double tau_w1minus = 70;
double tau_w2minus = 20;
double k_wminus = 65;
double u_wminus = 0.03;
double tau_wplus = 280;
double tau_fi = 0.11;
double tau_o1 = 6;
double tau_o2 = 6;
double tau_so1 = 43;
double tau_so2 = 0.2;
double k_so = 2;
double u_so = 0.65;
double tau_s1 = 2.7342;
double tau_s2 = 3;
double k_s = 2.0994;
double u_s = 0.9087;
double tau_si = 2.8723;
double tau_winf = 0.07;
double w_infstar = 0.94;
double D = 1;

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

double vnoinf(double x, double y)
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
void solve_tridiagonal_in_place_reusable(double * restrict const x, const size_t X, const double * restrict const a, const double * restrict const b, const double * restrict const c) {
    /* Allocate scratch space. */
    double * restrict const cprime = malloc(sizeof(double) * X);
    
    if (!cprime) {
        printf("Error: malloc failed.\n");
    }
    printf("a[0] as v = %f\n", a[0]);
    cprime[0] = c[0] / b[0];
    x[0] = x[0] / b[0];
    
    /* loop from 1 to X - 1 inclusive */
    for (size_t ix = 1; ix < X; ix++) {
        const double m = 1.0f / (b[ix] - a[ix] * cprime[ix - 1]);
        cprime[ix] = c[ix] * m;
        x[ix] = (x[ix] - a[ix] * x[ix - 1]) * m;
    }
    
    /* loop from X - 2 to 0 inclusive, safely testing loop end condition */
    for (size_t ix = X - 1; ix-- > 0; )
        x[ix] -= cprime[ix] * x[ix + 1];
    x[0] -= cprime[0] * x[1];
   printf("b[0] as x = %f\n", b[0]);
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

    // For a cell with 12000 time steps
    // Taxas de constantes no tempo
    /* double tauv, tauw, taus0, taus, tau0, vinf, winf, jfi, jso, jsi;
    double du_dt, dv_dt, ds_dt, dw_dt;

    // Variáveis
    double udef[12000], vdef[12000], wdef[12000], sdef[12000], tdef[12000];
    // Condicões iniciais
    double u = 0;
    double v = 1;
    double w = 1;
    double s = 0;
    double t = 0;

    udef[0] = 0;
    vdef[0] = v;
    wdef[0] = w;
    sdef[0] = s;
    tdef[0] = t; */
    /* int i;
    for (i = 1; i < 12000; i++)
    {
        if (i == 1000)
            u = 0.325;

        tauv = (1 - H(u, thetavlinha)) * tauv1 + H(u, thetavlinha) * tauv2;
        tauw = tauw1 + (tauw2 - tauw1) * (1 + tanh(kw * (u - uw))) / 2.0;
        taus0 = taus01 + (taus02 - taus01) * (1 + tanh(kso * (u - uso))) / 2.0;
        taus = (1 - H(u, thetaw)) * taus1 + H(u, thetaw) * taus2;
        tau0 = (1 - H(u, theta0)) * tau01 + H(u, theta0) * tau02;

        vinf = vnoinf(u, thetavlinha);
        winf = (1 - H(u, theta0)) * (1 - u / tauwinf) + H(u, theta0) * winfstar;

        jfi = -v * H(u, thetav) * (u - thetav) * (uu - u) / taufi;
        jso = (u - u0) * (1 - H(u, thetaw)) / tau0 + H(u, thetaw) / taus0;
        jsi = -H(u, thetaw) * w * s / tausi;

        du_dt = -(jfi + jso + jsi);
        dv_dt = (1 - H(u, thetav)) * (vinf - v) / tauv - H(u, thetav) * v / tauvmais;
        dw_dt = (1 - H(u, thetaw)) * (winf - w) / tauw - H(u, thetaw) * w / tauwmais;
        ds_dt = ((1 + tanh(ks * (u - us))) / 2 - s) / taus;

        u = u + delta_t * du_dt;
        v = v + delta_t * dv_dt;
        w = w + delta_t * dw_dt;
        s = s + delta_t * ds_dt;
        t = t + delta_t;

        udef[i] = u;
        vdef[i] = v;
        wdef[i] = w;
        sdef[i] = s;
        tdef[i] = t;
    } */
    // -----------------------------

    // Discretization
    int T = 200;
    double L = 1200;
    double delta_x = 1;
    double delta_t = 0.05;
    double delta_t_exp = 0.01;

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

    double u, v, w, s;
    double J_fi, J_so, J_si, J;
    double tau_vminus, tau_wminus, tau_so, tau_s, tau_o;
    double v_inf, w_inf;
    double du_dt, dv_dt, ds_dt, dw_dt;

    // Initial conditions
    u = 0;
    v = 1;
    w = 1;
    s = 0;

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
    a = malloc((N-1) * sizeof(double));     // subdiagonal
    b = malloc(N * sizeof(double));         // diagonal
    c = malloc((N-1) * sizeof(double));     // superdiagonal
    x = malloc(N * sizeof(double));


    for(int i = 0; i < N; i++)
    {
        if (i < N-1)
        {
            a[i] = -alpha;
            c[i] = -alpha;
        }
        a[i] = 1 + 2*alpha;
    }

    // Operator splitting method
    for (int n = 0; n < M - 1; n++)
    { 

        int n_exp, i;
        for (n_exp = 0; n_exp < M_exp; n_exp++)
        {

            for(i = 1; i < N-1; i++)
            {
                if(n == 1000)
                {
                    u = 0.325;
                }

                tau_vminus = (1-H(u, theta_vminus))*tau_v1minus + H(u, theta_vminus)*tau_v2minus;
                tau_wminus = tau_w1minus + (tau_w2minus-tau_w1minus)*(1+tanh(k_wminus*(u-u_wminus)))/2;
                tau_so = tau_so1 + (tau_so2-tau_so1)*(1+tanh(k_so*(u-u_so)))/2;
                tau_s = (1-H(u, theta_w))*tau_s1 + H(u, theta_w)*tau_s2;
                tau_o = (1-H(u, theta_o))*tau_o1 + H(u, theta_o)*tau_o2;

                J_fi = -v * H(u, theta_v) * (u-theta_v) * (u_u-u) / tau_fi;
                J_so = (u-u_o) * ((1-H(u, theta_w)) / tau_o) + (H(u, theta_w) / tau_so);
                J_si = -H(u, theta_w) * w * s / tau_si;
                J = J_fi + J_so + J_si;
                
                v_inf = vnoinf(u, theta_vminus);
                w_inf = (1-H(u, theta_o))*(1-u/tau_winf) + H(u, theta_o)*w_infstar;
                
                du_dt = - J;
                dv_dt = (1-H(u, theta_v))*(v_inf-v)/tau_vminus - H(u, theta_v)*v/tau_vplus;
                dw_dt = (1-H(u, theta_w))*(w_inf-w)/tau_wminus - H(u, theta_w)*w/tau_wplus;
                ds_dt = ((1+tanh(k_s*(u-u_s)))/2 - s)/tau_s;
                
                // Update variables
                u = u + du_dt * delta_t_exp;
                v = v + dv_dt * delta_t_exp;
                w = w + dw_dt * delta_t_exp;
                s = s + ds_dt * delta_t_exp;

            }



        }
        
        // REVER --------------------------------------------------

        for (i = 1; i < N - 1; i++)
        {
            if (n == 1000)
                u = 0.325; 

            int step_ni = n*N+i;
            int step_np1i = (n+1)*N+i;

            // Explicit part: J and ODEs (v, w, s)
            int i_exp;
            
            for(i_exp = 0; i_exp < M_exp; i_exp++)
            {
                if((n == 0 || n == 1) && i == 1)
                {
                    printf("\nAntes\n");
                    printf("u = %f\n", u);
                    printf("v = %f\n", v);
                    printf("w = %f\n", w);
                    printf("s = %f\n", s);
                }

                tau_vminus = (1-H(u, theta_vminus))*tau_v1minus + H(u, theta_vminus)*tau_v2minus;
                tau_wminus = tau_w1minus + (tau_w2minus-tau_w1minus)*(1+tanh(k_wminus*(u-u_wminus)))/2;
                tau_so = tau_so1 + (tau_so2-tau_so1)*(1+tanh(k_so*(u-u_so)))/2;
                tau_s = (1-H(u, theta_w))*tau_s1 + H(u, theta_w)*tau_s2;
                tau_o = (1-H(u, theta_o))*tau_o1 + H(u, theta_o)*tau_o2;

                J_fi = -v * H(u, theta_v) * (u-theta_v) * (u_u-u) / tau_fi;
                J_so = (u-u_o) * ((1-H(u, theta_w)) / tau_o) + (H(u, theta_w) / tau_so);
                J_si = -H(u, theta_w) * w * s / tau_si;
                J = J_fi + J_so + J_si;
                
                if((n == 0 || n == 1) && i == 1)
                {
                    printf("J_fi = %f; v = %f; H(u,vth) = %f; u = %f; tau_fi = %f  \n", J_fi, v, H(u, theta_v), u, tau_fi);
                    printf("J_so = %f; u = %f; H(u,thw) = %f; u_o = %f; tau_o = %f; tau_so = %f  \n", J_so, u, H(u, theta_w), u_o, tau_o, tau_so);
                    printf("J_si = %f; w = %f; s = %f; tau_si = %f\n", J_si, w, s, tau_si);
                    printf("J = %f\n", J);
                }
                
                v_inf = vnoinf(u, theta_vminus);
                w_inf = (1-H(u, theta_o))*(1-u/tau_winf) + H(u, theta_o)*w_infstar;
                
                du_dt = - J;
                dv_dt = (1-H(u, theta_v))*(v_inf-v)/tau_vminus - H(u, theta_v)*v/tau_vplus;
                dw_dt = (1-H(u, theta_w))*(w_inf-w)/tau_wminus - H(u, theta_w)*w/tau_wplus;
                ds_dt = ((1+tanh(k_s*(u-u_s)))/2 - s)/tau_s;

                if((n == 0 || n == 1) && i == 1)
                {
                    printf("dudt = %f\n", du_dt);
                }
                
                // Update variables
                u = u + du_dt * delta_t_exp;
                v = v + dv_dt * delta_t_exp;
                w = w + dw_dt * delta_t_exp;
                s = s + ds_dt * delta_t_exp;

                if((n == 0 || n == 1) && i == 1)
                {
                    printf("\nDepois\n");
                    printf("u = %f\n", u);
                    printf("v = %f\n", v);
                    printf("w = %f\n", w);
                    printf("s = %f\n", s);
                }
            }
            
            // Implicit part: PDE (u)
            x[0] = U[step_ni-1];
            x[1] = U[step_ni] + u;
            x[2] = U[step_ni+1];

            //printf("x[0] as v = %f\n", x[0]);

            solve_tridiagonal_in_place_reusable(x, 3, a, b, c);
            
            //if(n==1000)
                //printf("x[0] as x = %f\n", x[0]);

            U[step_np1i-1] = x[0];
            U[step_np1i] = x[1] + u;
            U[step_np1i+1] = x[2];

            //printf("u[stepnp1i-1] from x = %f\n", U[step_np1i-1]);
            
            u = U[step_np1i];
        }

        // Boundary Condition
        U[(n + 1) * N + 0] = U[(n + 1) * N + 1];
        U[(n + 1) * N + (N - 1)] = U[(n + 1) * N + (N - 2)];
    } 


    FILE *fp = NULL;
    fp = fopen("ceq-mm.dat", "w");
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
