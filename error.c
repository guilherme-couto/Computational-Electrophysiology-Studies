#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

int main()
{
    FILE *fp;
    fp = fopen("error.txt", "w");

    FILE *fp_1;
    fp_1 = fopen("mono_0.01.txt", "r");

    char *ptr;

    // Compare 0.01 and 0.02
    FILE *fp_2;
    fp_2 = fopen("mono_0.02.txt", "r");

    fprintf(fp, "Error (0.01 and 0.02):\n");

    double sum = 0;
    double sum_norm_01 = 0;
    double sum_abs = 0;

    for(int i = 0; i < 150*150; i ++)
    {
        char line_1[10];
        char line_2[10];

        fgets(line_1, 10, fp_1);
        fgets(line_2, 10, fp_2);

        double num_1 = strtod(line_1, &ptr);
        double num_2 = strtod(line_2, &ptr);

        double p = pow(num_2 - num_1, 2); 

        sum += p;
        sum_abs += fabs(num_2 - num_1);
        sum_norm_01 += pow(num_1, 2);
    }

    double rmse = sqrt(sum / (150*150));    // Root Mean Square Error
    double norm_01 = sqrt(sum_norm_01);    // Normalization
    double ed = sqrt(sum);    // Euclidean distance

    fprintf(fp, "Euclidean distance = %lf\n", ed);
    fprintf(fp, "RMSE = %lf\n", rmse);
    fprintf(fp, "Mean absolute error = %lf\n", sum_abs / (150*150));
    fprintf(fp, "ED / Norm_ref = %lf\n", (ed / norm_01));
    fprintf(fp, "Relative %% = %lf\n", (ed / norm_01) * 100);
    fprintf(fp, "\n\n");
    fclose(fp_2);
    fclose(fp_1);



    // Compare 0.01 and 0.03
    fp_1 = fopen("mono_0.01.txt", "r");
    FILE *fp_3;
    fp_3 = fopen("mono_0.03.txt", "r");

    fprintf(fp, "Error (0.01 and 0.03):\n");

    sum = 0;
    sum_abs = 0;

    for(int i = 0; i < 150*150; i ++)
    {
        char line_1[10];
        char line_3[10];

        fgets(line_1, 10, fp_1);
        fgets(line_3, 10, fp_3);

        double num_1 = strtod(line_1, &ptr);
        double num_3 = strtod(line_3, &ptr);

        double p = pow(num_3 - num_1, 2); 

        sum += p;
        sum_abs += fabs(num_3 - num_1);
    }

    rmse = sqrt(sum / (150*150));    // Root Mean Square Error
    ed = sqrt(sum);    // Euclidean distance

    fprintf(fp, "Euclidean distance = %lf\n", ed);
    fprintf(fp, "RMSE = %lf\n", rmse);
    fprintf(fp, "Mean absolute error = %lf\n", sum_abs / (150*150));
    fprintf(fp, "ED / Norm_ref = %lf\n", (ed / norm_01));
    fprintf(fp, "Relative %% = %lf\n", (ed / norm_01) * 100);
    fprintf(fp, "\n\n");
    fclose(fp_3);
    fclose(fp_1);



    // Compare 0.01 and 0.04
    fp_1 = fopen("mono_0.01.txt", "r");
    FILE *fp_4;
    fp_4 = fopen("mono_0.04.txt", "r");

    fprintf(fp, "Error (0.01 and 0.04):\n");

    sum = 0;
    sum_abs = 0;

    for(int i = 0; i < 150*150; i ++)
    {
        char line_1[10];
        char line_4[10];

        fgets(line_1, 10, fp_1);
        fgets(line_4, 10, fp_4);

        double num_1 = strtod(line_1, &ptr);
        double num_4 = strtod(line_4, &ptr);

        double p = pow(num_4 - num_1, 2); 

        sum += p;
        sum_abs += fabs(num_4 - num_1);
    }

    rmse = sqrt(sum / (150*150));    // Root Mean Square Error
    ed = sqrt(sum);    // Euclidean distance

    fprintf(fp, "Euclidean distance = %lf\n", ed);
    fprintf(fp, "RMSE = %lf\n", rmse);
    fprintf(fp, "Mean absolute error = %lf\n", sum_abs / (150*150));
    fprintf(fp, "ED / Norm_ref = %lf\n", (ed / norm_01));
    fprintf(fp, "Relative %% = %lf\n", (ed / norm_01) * 100);
    fprintf(fp, "\n\n");
    fclose(fp_4);
    fclose(fp_1);



    // Compare 0.01 and 0.05
    fp_1 = fopen("mono_0.01.txt", "r");
    FILE *fp_5;
    fp_5 = fopen("mono_0.05.txt", "r");

    fprintf(fp, "Error (0.01 and 0.05):\n");

    sum = 0;
    sum_abs = 0;

    for(int i = 0; i < 150*150; i ++)
    {
        char line_1[10];
        char line_5[10];

        fgets(line_1, 10, fp_1);
        fgets(line_5, 10, fp_5);

        double num_1 = strtod(line_1, &ptr);
        double num_5 = strtod(line_5, &ptr);

        double p = pow(num_5 - num_1, 2); 

        sum += p;
        sum_abs += fabs(num_5 - num_1);
    }

    rmse = sqrt(sum / (150*150));    // Root Mean Square Error
    ed = sqrt(sum);    // Euclidean distance

    fprintf(fp, "Euclidean distance = %lf\n", ed);
    fprintf(fp, "RMSE = %lf\n", rmse);
    fprintf(fp, "Mean absolute error = %lf\n", sum_abs / (150*150));
    fprintf(fp, "ED / Norm_ref = %lf\n", (ed / norm_01));
    fprintf(fp, "Relative %% = %lf\n", (ed / norm_01) * 100);
    fprintf(fp, "\n\n");
    fclose(fp_5);

    fclose(fp_1);

    // Compare 0.01 (adi) and 0.01 (exp)
    fp_1 = fopen("mono_0.01.txt", "r");
    FILE *fp_exp;
    fp_exp = fopen("mono_exp_0.01.txt", "r");

    fprintf(fp, "Error (0.01 adi and 0.01 explicit):\n");

    sum = 0;
    sum_abs = 0;

    for(int i = 0; i < 150*150; i ++)
    {
        char line_1[10];
        char line_exp[10];

        fgets(line_1, 10, fp_1);
        fgets(line_exp, 10, fp_5);

        double num_1 = strtod(line_1, &ptr);
        double num_exp = strtod(line_exp, &ptr);

        double p = pow(num_exp - num_1, 2); 

        sum += p;
        sum_abs += fabs(num_exp - num_1);
    }

    rmse = sqrt(sum / (150*150));    // Root Mean Square Error
    ed = sqrt(sum);    // Euclidean distance

    fprintf(fp, "Euclidean distance = %lf\n", ed);
    fprintf(fp, "RMSE = %lf\n", rmse);
    fprintf(fp, "Mean absolute error = %lf\n", sum_abs / (150*150));
    fprintf(fp, "ED / Norm_ref = %lf\n", (ed / norm_01));
    fprintf(fp, "Relative %% = %lf\n", (ed / norm_01) * 100);
    fprintf(fp, "\n\n");
    fclose(fp_exp);

    fclose(fp_1);

    fclose(fp);

    return 0;
}