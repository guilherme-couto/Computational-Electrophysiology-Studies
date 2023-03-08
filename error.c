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
    char *ptr;

    // Compare EXP with dt = 0.01 and ADI with dt = 0.01
    fp_1 = fopen("mono_exp_0.01.txt", "r");
    FILE *fp_adi;
    fp_adi = fopen("mono_0.01.txt", "r");

    fprintf(fp, "Error (EXP with dt = 0.01 and ADI with dt = 0.01):\n");

    double sum = 0;
    double sum_norm_01 = 0;
    double sum_abs = 0;

    for(int i = 0; i < 150*150; i ++)
    {
        char line_1[10];
        char line_adi[10];
        fgets(line_1, 10, fp_1);
        fgets(line_adi, 10, fp_adi);
        double num_1 = strtod(line_1, &ptr);
        double num_exp = strtod(line_adi, &ptr);

        double p = pow(num_exp - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_exp - num_1);
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
    
    fclose(fp_adi);
    fclose(fp_1);


    // Compare EXP with dt = 0.01 and ADI dt = 0.02
    fp_1 = fopen("mono_exp_0.01.txt", "r");
    FILE *fp_2;
    fp_2 = fopen("mono_0.02.txt", "r");

    fprintf(fp, "Error (EXP with dt = 0.01 and ADI dt = 0.02):\n");

    sum = 0;
    sum_abs = 0;

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
    }

    rmse = sqrt(sum / (150*150));    // Root Mean Square Error
    ed = sqrt(sum);    // Euclidean distance

    fprintf(fp, "Euclidean distance = %lf\n", ed);
    fprintf(fp, "RMSE = %lf\n", rmse);
    fprintf(fp, "Mean absolute error = %lf\n", sum_abs / (150*150));
    fprintf(fp, "ED / Norm_ref = %lf\n", (ed / norm_01));
    fprintf(fp, "Relative %% = %lf\n", (ed / norm_01) * 100);
    fprintf(fp, "\n\n");
    
    fclose(fp_2);
    fclose(fp_1);


    // Compare EXP with dt = 0.01 and ADI dt = 0.03
    fp_1 = fopen("mono_exp_0.01.txt", "r");
    FILE *fp_3;
    fp_3 = fopen("mono_0.03.txt", "r");

    fprintf(fp, "Error (EXP with dt = 0.01 and ADI dt = 0.03):\n");

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


    // Compare EXP with dt = 0.01 and ADI dt = 0.04
    fp_1 = fopen("mono_exp_0.01.txt", "r");
    FILE *fp_4;
    fp_4 = fopen("mono_0.04.txt", "r");

    fprintf(fp, "Error (EXP with dt = 0.01 and ADI dt = 0.04):\n");

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


    // Compare EXP with dt = 0.01 and ADI dt = 0.05
    fp_1 = fopen("mono_exp_0.01.txt", "r");
    FILE *fp_5;
    fp_5 = fopen("mono_0.05.txt", "r");

    fprintf(fp, "Error (EXP with dt = 0.01 and ADI dt = 0.05):\n");

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


    // Compare EXP with dt = 0.01 and ADI dt = 0.1
    fp_1 = fopen("mono_exp_0.01.txt", "r");
    FILE *fp_10;
    fp_10 = fopen("mono_0.1.txt", "r");

    fprintf(fp, "Error (EXP with dt = 0.01 and ADI dt = 0.1):\n");

    sum = 0;
    sum_abs = 0;

    for(int i = 0; i < 150*150; i ++)
    {
        char line_1[10];
        char line_10[10];
        fgets(line_1, 10, fp_1);
        fgets(line_10, 10, fp_10);
        double num_1 = strtod(line_1, &ptr);
        double num_10 = strtod(line_10, &ptr);

        double p = pow(num_10 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_10 - num_1);
    }

    rmse = sqrt(sum / (150*150));    // Root Mean Square Error
    ed = sqrt(sum);    // Euclidean distance

    fprintf(fp, "Euclidean distance = %lf\n", ed);
    fprintf(fp, "RMSE = %lf\n", rmse);
    fprintf(fp, "Mean absolute error = %lf\n", sum_abs / (150*150));
    fprintf(fp, "ED / Norm_ref = %lf\n", (ed / norm_01));
    fprintf(fp, "Relative %% = %lf\n", (ed / norm_01) * 100);
    fprintf(fp, "\n\n");
    
    fclose(fp_10);
    fclose(fp_1);


    // Compare EXP with dt = 0.01 and ADI dt = 0.2
    fp_1 = fopen("mono_exp_0.01.txt", "r");
    FILE *fp_20;
    fp_20 = fopen("mono_0.2.txt", "r");

    fprintf(fp, "Error (EXP with dt = 0.01 and ADI dt = 0.2):\n");

    sum = 0;
    sum_abs = 0;

    for(int i = 0; i < 150*150; i ++)
    {
        char line_1[10];
        char line_20[10];
        fgets(line_1, 10, fp_1);
        fgets(line_20, 10, fp_20);
        double num_1 = strtod(line_1, &ptr);
        double num_20 = strtod(line_20, &ptr);

        double p = pow(num_20 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_20 - num_1);
    }

    rmse = sqrt(sum / (150*150));    // Root Mean Square Error
    ed = sqrt(sum);    // Euclidean distance

    fprintf(fp, "Euclidean distance = %lf\n", ed);
    fprintf(fp, "RMSE = %lf\n", rmse);
    fprintf(fp, "Mean absolute error = %lf\n", sum_abs / (150*150));
    fprintf(fp, "ED / Norm_ref = %lf\n", (ed / norm_01));
    fprintf(fp, "Relative %% = %lf\n", (ed / norm_01) * 100);
    fprintf(fp, "\n\n");
    
    fclose(fp_20);
    fclose(fp_1);


    // Compare EXP with dt = 0.01 and ADI dt = 0.5
    fp_1 = fopen("mono_exp_0.01.txt", "r");
    FILE *fp_50;
    fp_50 = fopen("mono_0.5.txt", "r");

    fprintf(fp, "Error (EXP with dt = 0.01 and ADI dt = 0.5):\n");

    sum = 0;
    sum_abs = 0;

    for(int i = 0; i < 150*150; i ++)
    {
        char line_1[10];
        char line_50[10];
        fgets(line_1, 10, fp_1);
        fgets(line_50, 10, fp_50);
        double num_1 = strtod(line_1, &ptr);
        double num_50 = strtod(line_50, &ptr);

        double p = pow(num_50 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_50 - num_1);
    }

    rmse = sqrt(sum / (150*150));    // Root Mean Square Error
    ed = sqrt(sum);    // Euclidean distance

    fprintf(fp, "Euclidean distance = %lf\n", ed);
    fprintf(fp, "RMSE = %lf\n", rmse);
    fprintf(fp, "Mean absolute error = %lf\n", sum_abs / (150*150));
    fprintf(fp, "ED / Norm_ref = %lf\n", (ed / norm_01));
    fprintf(fp, "Relative %% = %lf\n", (ed / norm_01) * 100);
    fprintf(fp, "\n\n");
    
    fclose(fp_50);
    fclose(fp_1);

    fclose(fp);
    return 0;
}