#include "opti_method.h"

#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>

int main(){
    double lambda0, lambda1, x0[2], x1[2], xHat[2], tol, lamb_opt;
    int maxIter;
    lambda0 = 0.0;
    lambda1 = 1.0;

    x0[0] = 0;
    x0[1] = 1;
    x1[0] = 1;
    x1[1] = 0;
    xHat[0] = 1;
    xHat[1] = 1;
    tol = 0.001;
    maxIter = 10;
    lamb_opt = secant_2D(lambda0, lambda1, x0, x1, xHat, tol, maxIter);
    printf("\nx0: %f, %f", x0[0], x0[1]);
    printf("\nx1: %f, %f", x1[0], x1[1]);
    printf("Lamba 1: %f \n", lambda1);
    printf("Lamba 2: %f \n", lambda0);
    printf("Optimum lambda: %f \n", lamb_opt);

    x0[0] = 1;
    x0[1] = 2;
    x1[0] = 2;
    x1[1] = 1;
    xHat[0] = 1;
    xHat[1] = 1;
    tol = 0.001;
    maxIter = 10;
    lamb_opt = secant_2D(lambda0, lambda1, x0, x1, xHat, tol, maxIter);
    printf("\nx0: %f, %f", x0[0], x0[1]);
    printf("\nx1: %f, %f", x1[0], x1[1]);
    printf("Lamba 1: %f \n", lambda1);
    printf("Lamba 2: %f \n", lambda0);
    printf("Optimum lambda: %f \n", lamb_opt);


    x0[0] = 1;
    x0[1] = 0;
    x1[0] = 0;
    x1[1] = 1;
    xHat[0] = 1;
    xHat[1] = 1;
    tol = 0.0000001;
    maxIter = 10;
    lamb_opt = secant_2D(lambda0, lambda1, x0, x1, xHat, tol, maxIter);
    printf("\nx0: %f, %f", x0[0], x0[1]);
    printf("\nx1: %f, %f", x1[0], x1[1]);
    printf("Lamba 1: %f \n", lambda1);
    printf("Lamba 2: %f \n", lambda0);
    printf("Optimum lambda: %f \n", lamb_opt);

    // To test that it has the right convergence rate
    int N = 25;
    double *conv;
    conv = malloc(N*sizeof(double));

    for(int i=0; i<N; i++){
        maxIter = 1 + i;
        conv[i] = (double)secant_2D(lambda0, lambda1, x0, x1, xHat, tol, maxIter);
        printf("%f \n", conv[i]);
    }

    FILE *fp;
    fp = fopen("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/Outputs/convSec.bin", "wb");
    fwrite(conv, sizeof(double), N, fp);
    fclose(fp);

    free(conv);

}


