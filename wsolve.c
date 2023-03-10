#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include "wsolve.h"

#define MAX_N       128
#define CONFIGFILE  "wsolve.conf"

// Globals
const char def_configfile[] = CONFIGFILE;
const char help_text[] = "wsolve [options] <infile> <outfile>\n"\
"\n"\
"The WSOLVE binary parses raw data from a Spinning Disc Langmuir Probe\n"\
"to establish spatially resolved wire current density, Ibar(x,y).  The\n"\
"output of WSOLVE is a set of complex-valued coefficients for a Fourier\n"\
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
"<infile>\n"\
"  Raw data are read in double-precision floating point quartets from a\n"\
"  data file. A quartet includes the wire radius, R, the X and Y location\n"\
"  of the disc center in the domain, the disc angle in radians, THETA,\n"\
"  and the measured wire current in that configuration, I. The R,X,Y,THETA,I"\
"  groups repeat in the file with no header, footer, and with no separation,\n"\
"      ...\n"\
"    [R     double]\n"\
"    [X     double]\n"\
"    [Y     double]\n"\
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
"    [I0_real    double]\n"\
"    [I0_imag    double]\n"\
"\n"\
"  Each real and imaginary coefficient pair corresponds to an m- and n-index\n"\
"  in the expansion, above.  The indices specify the wavenumber in the x- and\n"\
"  y-axes.  The last value, I0, is an offset current that is found present at\n"\
"  all wire locations. This is not attributed to plasma currents, but to errors\n"\
"  in the wire current measurement calibration.\n"\
"\n"\
"  Because m varies from -Nx to +Nx, n varies from -Ny to +Ny, and the total\n"\
"  solution includes all possible combinations of these, there are (2Nx+1)(2Ny+1)\n"\
"  complex-valued coefficients.  Each is stored with its real and imaginary \n"\
"  parts as double-precision floats in that order.\n"\
"\n"\
"  The coefficients are stored in order with m increasing first and then n. So\n"\
"  if Nx=Ny=1, the nine coefficients would appear in order (c_mn):\n"\
"    c_-1-1 c0-1 c1-1 c-10 c00 c10 c-11 c01 c11\n"\
"\n"\
"*** CONFIGURATION ***\n"\
"By default, the configuration file is \"wsolve.conf\" in the working dir.\n"\
"The configuration file is a plain text file that is used to define the\n"\
"shape of the computational domain, the number of wavenumbers to use in\n"\
"the model, and the number of threads to spawn for the job. It should be\n"\
"in the working directory. Its contents must appear:\n"\
"\n"
"    nthread <int>\n"\
"    Nx <int>\n"\
"    Ny <int>\n"\
"    Lx <float>\n"\
"    Ly <float>\n"\
"    xshift <float>\n"\
"    yshift <float>\n"\
"\n"\
"The \"nthread\" integer specifies the number of threads that should be\n"\
"spawned to perform the computation.  The Nx and Ny integers specify the\n"\
"numebr of wavenumbers on each axis.  Lx and Ly specify the rectangular\n"\
"domain width and height.\n"\
"\n"\
"Finally, xshift and yshift specify an offset to add to all x- and y-\n"\
"coordinates from the data files. This has the effect of shifting the\n"\
"rectangular domain relative to the data.\n"\
"\n"\
"To simplify the code, the configuration file format is extremely strict.\n"\
"All characters are case sensitive, and the order must not be changed.\n"\
"However, the amount of whitespace may be changed and any text beyond the\n"\
"dshift parameters will not be read.\n"\
"\n"\
"*** OPTIONS ***\n"\
"-c <config_file>\n"\
"  Override the default configuration file: wsolve.conf\n"\
"\n"\
"-h\n"\
"  Print this help and exit.\n"\
"\n"\
"-q\n"\
"  Run quietly; disables printing status messages to stdout.\n"\
"\n"\
"(c)2023 Christopher R. Martin\n"\
"\n";



int main(int argc, char *argv[]){
    unsigned int Nx,Ny,nthread;
    double Lx,Ly;
    double xshift, yshift;
    char verbose_f;
    char *configfile, *infile, *outfile;
    char ch;
    int ii, jj, kk;
    WireSlice_t ws;
    int err;
    FILE * fd;
    time_t start, now;

    // Set default settings
    configfile = NULL;
    verbose_f = 1;

    // Parse the options
    while((ch = getopt(argc, argv, "hqc:")) >= 0){
        switch(ch){
        case 'h':
            printf(help_text);
            return 0;
        case 'q':
            verbose_f = 0;
        break;
        case 'c':
            configfile = optarg;
        break;
        default:
            fprintf(stderr, "WSOLVE: Unrecognized option %c\n", (char) ch);
            return -1;
        break;
        }
    }

    // Find the infile and outfile arguments
    // If there aren't two arguments left, raise an error
    if(argc - optind != 2){
        fprintf(stderr, "WSOLVE: Requires two non-option arguments - type \"wsolve -h\" for more information\n");
        return -1;
    }
    infile = argv[optind];
    outfile = argv[optind+1];


    // Read the config file.
    // Open the configuration file
    if(!configfile)
        configfile = def_configfile;
    fd = fopen(configfile, "r");
    if(!fd){
        fprintf(stderr, "WSOLVE: Failed to open configuration file:\n  %s\n", configfile);
        return -1;
    }
    
    // Set the defaults
    Nx = 0;
    Ny = 0;
    Lx = 0;
    Ly = 0;
    xshift = 0;
    yshift = 0;
    nthread = 1;
    // Read the file
    err = fscanf(fd, "nthread %d Nx %d Ny %d Lx %lf Ly %lf xshift %lf yshift %lf",
            &nthread, &Nx, &Ny, &Lx, &Ly, &xshift, &yshift);
    fclose(fd);
    fd = NULL;
    // Check that all the parameters were found
    if(err != 7){
        fprintf(stderr, "Configuration syntax error in %s. Call \"wsolve -h\" for help.\n", configfile);
        return -1;
    }
    err = 0;
    // Test the parameters for legal values
    if(Nx > MAX_N || Nx == 0){
        fprintf(stderr, "WSOLVE.CONF value error. Nx is illegal: %d\n", Nx);
        return -1;
    }else if(Ny > MAX_N || Ny == 0){
        fprintf(stderr, "WSOLVE.CONF value error. Ny is illegal: %d\n", Ny);
        return -1;
    }else if(Lx <= 0){
        fprintf(stderr, "WSOLVE.CONF value error. Lx is illegal: %lf\n", Lx);
        return -1;
    }else if(Ly <= 0){
        fprintf(stderr, "WSOLVE.CONF value error. Ly is illegal: %lf\n", Ly);
        return -1;
    }
    
        
    err = ws_init(&ws, Nx, Ny, Lx, Ly, verbose_f);
    err = err || ws_shift(&ws, xshift, yshift);
    err = err || ws_read(&ws, infile, nthread);
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
    err = err || ws_write(&ws, outfile);

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



