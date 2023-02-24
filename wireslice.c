#include <stdio.h>
#include <unistd.h>
#include "wireslice.h"

#define MAX_FILENAME 128
#define CONFIGFILE "wireslice.conf"

// Globals
const char [] help_text = "wireslice <datafile>\n"\
"\n"\
"The WIRESLICE binary parses raw data from a spinning disc Langmuir probe\n"\
"to establish spatially resolved wire current density.\n"\
"\n"\
"<datafile>\n"\
"  Raw data are read in double-precision floating point quartets from a\n"\
"  data file. A quartet includes the wire radius, R, the disc distance\n"\
"  from the edge of the measurement domain, D, the disc angle in radians,\n"\
"  THETA, and the measured wire current in that configuration, I. The \n"\
"  R,D,THETA,I quartets repeat in the file with no header, footer, and\n"\
"  with no separation,\n"\
"    ...[R][D][THETA][I][R][D][THETA][I]...\n"\
"\n"\
"wireslice.conf\n"\
"  The configuration file is a plain text file that is used to define the\n"\
"  shape of the computational domain, the number of wavenumbers to use in\n"\
"  the model, and the number of threads to spawn for the job. It should be\n"\
"  in the working directory."\
"\n"\


int main(int argc, char *argv[]){
    unsigned int Nx,Ny;
    double Lx,Ly;
    WireSlice_t ws;
    int err;
    FILE * fd;

    // Read in the configuration file
    if(!(fd = fopen(CONFIGFILE, "r"))){
        fprintf(stderr, "WIRESLICE: Failed to open configuration file:\n  %s\n", CONFIGFILE);
        return -1;
    }
    

    // Verify that the call signature appears correct
    if(argc != 2){
        printf(help_text);
        return -1;
    }
    
    err = ws_init(&ws, 20, 20, 1, 1);
    err = err || ws_read(&ws, "test.wf", 8);
    err = err || ws_solve(&ws);
    err = err || ws_write(&ws, "output.wf");
    ws_destruct(&ws);
    return err;
}



