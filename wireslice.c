#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include "wireslice.h"

#define MAX_N       128
#define CONFIGFILE  "wireslice.conf"
#define LOGFILE     "wireslice.log"
#define OUTFILE     "wireslice.out"

// Globals
const char help_text[] = "wireslice <infile> <outfile>\n"\
"\n"\
"The WIRESLICE binary parses raw data from a Spinning Disc Langmuir Probe\n"\
"to establish spatially resolved wire current density, Ibar(x,y).  The\n"\
"output of WIRESLICE is a set of complex-valued coefficients for a Fourier\n"\
"series on x and y that optimally match the raw measurements.\n"\
"\n"\
"            Nx     Ny              /       ->   -> \\ \n"\
"    Ibar = sum    sum    c_m,n exp | 2pi j nu . x  | \n"\
"           m=-Nx  n=-Ny            \\               / \n"\
"\n"\
"    ->    m  ^      n  ^\n"\
"    nu = --- i  +  --- j\n"\
"          Lx        Ly\n"\
"\n"\
"Here, Lx and Ly are the horizontal and vertical size of the rectangular\n"\
"domain.  The indices, m and n, form a wave number vector, nu.  The\n"\
"complexity that can be represented by the expansion is determined by Nx and\n"\
"Ny.\n"\
"\n"\
"wireslice.conf\n"\
"  The configuration file is a plain text file that is used to define the\n"\
"  shape of the computational domain, the number of wavenumbers to use in\n"\
"  the model, and the number of threads to spawn for the job. It should be\n"\
"  in the working directory. Its contents must appear:\n"\
"\n"
"    nthread <int>\n"\
"    x <int> <float>\n"\
"    y <int> <float>\n"\
"\n"\
"  The \"nthread\" integer specifies the number of threads that should be\n"\
"  spawned to perform the computation.  The integers specify the numebr of\n"\
"  wavenumbers on each axis.\n"\
"\n"\
"  To simplify the code, the configuration file format is extremely strict.\n"\
"  All characters should be lower case, and the order must not be changed.\n"\
"  However, the amount of whitespace may be changed and any text beyond the\n"\
"  y parameters will not be read.\n"\
"\n"\
"<infile>\n"\
"  Raw data are read in double-precision floating point quartets from a\n"\
"  data file. A quartet includes the wire radius, R, the disc distance\n"\
"  from the edge of the measurement domain, D, the disc angle in radians,\n"\
"  THETA, and the measured wire current in that configuration, I. The \n"\
"  R,D,THETA,I quartets repeat in the file with no header, footer, and\n"\
"  with no separation,\n"\
"      ...\n"\
"    [R     double]\n"\
"    [D     double]\n"\
"    [THETA double]\n"\
"    [I     double]\n"\
"      ...\n"\
"\n"\
"<outfile>\n"\
"  The results of the transform is a binary file containing the domain\n"\
"  size, number of coefficients, and real and imaginary values for each \n"\
"  coefficient.  The header is two integers and two doubles describing\n"\
"  the domain:\n"\
"    [Nx uint32]\n"\
"    [Ny uint32]\n"\
"    [Lx double]\n"\
"    [Ly double]\n"\
"      ...\n"\
"    [c_k-1_real double]\n"\
"    [c_k-1_imag double]\n"\
"    [c_k_real   double]\n"\
"    [c_k_imag   double]\n"\
"    [c_k+1_real double]\n"\
"    [c_k+1_imag double]\n"\
"      ...\n"\
"\n"\
"  Each real and imaginary coefficient pair corresponds to an m- and n-index\n"\
"  in the expansion, above.  The indices specify the wavenumber in the x- and\n"\
"  y-axes.\n"\
"\n"\
"  Because m varies from -Nx to +Nx, n varies from -Ny to +Ny, and the total\n"\
"  solution includes all possible combinations of these, there are (2Nx+1)(2Ny+1)\n"\
"  complex-valued coefficients.  Each is stored with its real and imaginary \n"\
"  parts as double-precision floats in that order.\n"\
"\n"\
"  The coefficients are stored in order with m increasing first and then n. So\n"\
"  if Nx=Ny=1, the nine coefficients would appear in order (c_mn):\n"\
"    c_-1-1 c0-1 c1-1 c-10 c00 c10 c-11 c01 c11\n"\
"\n";


int main(int argc, char *argv[]){
    unsigned int Nx,Ny,nthread;
    double Lx,Ly;
    int ii, jj, kk;
    WireSlice_t ws;
    int err;
    FILE * fd;
    time_t start, now;

    // Verify that the call signature appears correct
    if(argc != 3){
        printf(help_text);
        return -1;
    }

    // Open the configuration file
    if(!(fd = fopen(CONFIGFILE, "r"))){
        fprintf(stderr, "WIRESLICE: Failed to open configuration file:\n  %s\n", CONFIGFILE);
        return -1;
    }
    
    // Set the defaults
    Nx = 0;
    Ny = 0;
    Lx = 0;
    Ly = 0;
    nthread = 1;
    // Read the file
    err = fscanf(fd, "nthread %d x %d %lf y %d %lf",
            &nthread, &Nx, &Lx, &Ny, &Ly);
    fclose(fd);
    // Check that all the parameters were found
    if(err != 5){
        fprintf(stderr, "WIRESLICE.CONF syntax error. Call wireslice with no arguments for help.\n");
        return -1;
    }
    // Test the parameters for legal values
    if(Nx > MAX_N || Nx == 0){
        fprintf(stderr, "WIRESLICE.CONF value error. Nx is illegal: %d\n", Nx);
        return -1;
    }else if(Ny > MAX_N || Ny == 0){
        fprintf(stderr, "WIRESLICE.CONF value error. Ny is illegal: %d\n", Ny);
        return -1;
    }else if(Lx <= 0){
        fprintf(stderr, "WIRESLICE.CONF value error. Lx is illegal: %lf\n", Lx);
        return -1;
    }else if(Ly <= 0){
        fprintf(stderr, "WIRESLICE.CONF value error. Ly is illegal: %lf\n", Ly);
        return -1;
    }
    
    start = time(NULL);
    printf("[%6.1fs] Initializing\n", difftime(time(NULL), start));
        
    err = ws_init(&ws, Nx, Ny, Lx, Ly);
    
    printf("[%6.1fs] WS_INIT RC: %d\n", 
            difftime(time(NULL), start),
            err);
    // WS_READ prints its own status messages to stdout
    err = err || ws_read(&ws, argv[1], nthread);

    printf("[%6.1fs] WS_READ RC: %d\n",
            difftime(time(NULL), start),
            err);
    /*
     * This code prints the A matrix and B vectors for debugging
     * 
    for(ii=0;ii<ws.ncoef;ii++){
        for(jj=0;jj<ws.ncoef;jj++){
            if(jj >= ii){
                kk = ii + jj*(jj+1)/2;
                printf("%7.1lf+%7.1lfj, ", creal(ws.AP[kk]), cimag(ws.AP[kk]));
            }else{
                kk = jj + ii*(ii+1)/2;
                printf("%7.1lf+%7.1lfj, ", creal(ws.AP[kk]), -cimag(ws.AP[kk]));
            }
        }
        printf("\n");
    }
    
    for(ii=0;ii<ws.ncoef;ii++){
        printf("%lf + %lfj\n", creal(ws.B[ii]), cimag(ws.B[ii]));
    }
    */
    
    err = err || ws_solve(&ws);

    printf("[%6.1fs] WS_SOLVE RC: %d\n",
            difftime(time(NULL), start),
            err);
                
    err = err || ws_write(&ws, argv[2]);

    printf("[%6.1fs] WS_WRITE RC: %d\n  Done.\n",
            difftime(time(NULL), start),
            err);

    /*
     * This code prints the solution vector for debugging
     * 
    for(ii=0;ii<ws.ncoef;ii++){
        printf("%lf + %lfj\n", creal(ws.B[ii]), cimag(ws.B[ii]));
    }
    */
    
    ws_destruct(&ws);
    
    return err;
}



