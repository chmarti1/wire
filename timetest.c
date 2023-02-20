#include <complex.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

#define REPEAT 10000
#define DIVS 1000
#define DTHETA 2*M_PI / DIVS

double complex ejth(double theta){
    double c;
    c = cos(theta);
    return CMPLX(c, sqrt(1-c*c));
}

int main(int argc, char *argv[]){
    int ii;
    double theta;
    double complex ee;
    clock_t tt;
    
    tt = clock();
    for(ii=0;ii<REPEAT;ii++){
        for(theta = 0; theta<2*M_PI; theta+=DTHETA){
            ee = ejth(theta);
        }
    }
    tt = clock() - tt;
    printf("EJTH: %lf\n", ((double) tt) / REPEAT / DIVS);
    
    tt = clock();
    for(ii=0;ii<REPEAT;ii++){
        for(theta = 0; theta<2*M_PI; theta+=DTHETA){
            ee = cexp(I * theta);
        }
    }
    tt = clock() - tt;
    printf("CEXP: %lf\n", ((double) tt) / REPEAT / DIVS);
    return 0;
}
