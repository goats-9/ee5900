#include <stdio.h>
#include <stdlib.h>
#define N (int)1e5

/**
 * Quantize a number x given precision prec
 */
double quantize(double x, int prec) { 
    if (x != -1) {
        long long X = x*(1LL<<prec);
        double xq = (X*1.0)/(1LL<<prec);
        return xq;
    }
    return x;
}

/**
 * Direct Form II Filter
 */
double *ap_df2_filter(double a, int prec) { 
    // Create input and output arrays
    double x[N] = {0};
    x[0] = 1;
    double v[N] = {0};
    v[0] = 1;
    double *y = (double *) malloc(N*sizeof(double));
    double A = quantize(a, prec);
    double Am = quantize(-a, prec);
    y[0] = Am;
    // Accumulate output
    for (int i = 1; i < N; i++) { 
        v[i] = x[i] + quantize(A*v[i-1], prec);
        y[i] = quantize(Am*v[i], prec) + v[i-1];
    }
    return y;
}

/**
 * Reduced Form Filter
 */
double *ap_red_filter(double a, int prec) { 
    // Create input and output arrays
    double x[N] = {0};
    x[0] = 1;
    double *y = (double *) malloc(N*sizeof(double));
    double A = quantize(a, prec);
    double Am = quantize(-a, prec);
    y[0] = Am;
    // Accumulate output
    for (int i = 1; i < N; i++) { 
        y[i] = quantize(A*(y[i-1]-x[i]), prec) + x[i-1];
    }
    return y;
}

int main() {
    double a = 0.1;
    // Quantized outputs for Direct Form II filter
    double *yd1 = ap_df2_filter(a, -1), *yd2 = ap_df2_filter(a, 2);
    // Quantized outputs for Reduced filter
    double *yr1 = ap_red_filter(a, -1), *yr2 = ap_red_filter(a, 2);
    // Compute differences
    double ed[N], er[N];
    for (int i = 0; i < N; i++) { 
        ed[i] = yd2[i] - yd1[i];
        er[i] = yr2[i] - yr1[i];
    }
    // Write outputs to a file
    FILE *dfile = fopen("df2_error.txt", "w");
    FILE *rfile = fopen("red_error.txt", "w");
    for (int i = 0; i < N; i++) { 
        fprintf(dfile, "%lf\n", ed[i]);
        fprintf(rfile, "%lf\n", er[i]);
    }
    fclose(rfile);
    fclose(dfile);
    return 0;
}
