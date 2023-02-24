#include <complex.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

#define REPEAT 10000
#define DIVS 1000
#define DTHETA 2*M_PI / DIVS

// calculate exp(j*theta) efficiently
double complex ejtheta(double theta){
    double c_th, mod_theta;
    c_th = cos(theta);
    mod_theta = theta / (2*M_PI);
    mod_theta -= floor(mod_theta);
    if(mod_theta >= 0.5)
        return CMPLX(c_th, -sqrt(1. - c_th*c_th));
    return CMPLX(c_th, sqrt(1. - c_th*c_th));
}

int main(int argc, char *argv[]){
    int ii;
    double theta;
    double complex ee;
    clock_t tt;
    
    for(theta = -4*M_PI; theta<-2*M_PI; theta+=10*DTHETA){
        ee = ejtheta(theta);
        printf("%6.4lf : %7.4lf+%7.4lfj\n", theta, creal(ee), cimag(ee));
    }    
    
    tt = clock();
    for(ii=0;ii<REPEAT;ii++){
        for(theta = 0; theta<2*M_PI; theta+=DTHETA){
            ee = ejtheta(theta);
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
