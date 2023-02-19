#include <complex.h>
#include <lapacke.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>

// XXX - may need the -fopenmp compiler directive for parallel processing? 

#define False 0                 // Boolean constants
#define True 1                  //
#define Element double complex  // Coefficient and matrix data type
#define BUFFER_SIZE 128        // Lines to buffer from the source data
#define MAX_THREADS 12          // Maximum worker threads to spawn in ws_read()

/* The WireSlice struct
 *  Wire slice structs contain the complex-valued matrix problem for
 * expressing
 * 
 * Sets up a linear system to solve for C coefficients in the expansion
 *               Nx     Ny               /       ->  ->    \
 *      Ibar =  sum    sum      c    exp | 2j pi X . NU    |
 *             m=-Nx  n=-Ny      m,n     \             m,n /
 * 
 * where
 *      ->  ->          m*x       n*y
 *      X . NU      =  -----  +  -----
 *            m,n        Lx        Ly
 * 
 * given a series of measurements of Ibar integrated along line segments
 * passing through the domain.
 * 
 * The problem is expressed as a complex-valued system of linear equa-
 * tions,
 * 
 *      B = AP * C
 * 
 * The A-matrix is Hermitian (A = A*), so we have borrowed the "AP" 
 * notation used by LAPACK to emphasize that only the upper triangle is 
 * stored in memory.
 * 
 * WS_OPEN(), WS_READ(), and WS_CLOSE()
 * manage the file that contains the raw data that construct the problem.
 * The binary file should contain groups of four double-precision values
 * in order: (R, D, THETA, I).  R is the wire radius, D is the disc 
 * distance from the y-axis, THETA is the wire angle, and I is the 
 * measured wire current. 
 * 
 * WS_READ() is designed to be called multiple times in parallel by 
 * 
 */

typedef struct WireSlice {
    Element * AP;
    Element * B;
    Element * C;
    unsigned int Nx, Ny;
    double Lx, Ly;
    unsigned int ncoef;
    unsigned int nAP;
    unsigned int ndata;
    unsigned char _halt;
    unsigned char _alock;
    unsigned char _flock;
    FILE * target;
} WireSlice_t;

/* Prototypes
 *  The functions that operate on the WireSlice are:
 *  ws_init     : Initialize the WireSlice struct with wavenumber and dimensions
 *  ws_destruct : Free allocated memory and close the file (if open)
 *  ws_read     : Read all data into the struct
 *  ws_solve    : Solve the problem defined by ws_init and ws_read
 *  ws_write    : Write the solution results to a file
 * 
 */
 
/* WS_INIT  -   Initialize the WireSlice struct with wavenumber and dimensions
 *  ws      :   WireSlice struct
 *  Nx, Ny  :   Number of x- and y-wavenumbers
 *  Lx, Ly  :   Domain length in x- and y-axes
 * 
 * Returns 0 on success
 * Returns -1 on memory allocation failure
 */
int ws_init(WireSlice_t * ws, unsigned int Nx, unsigned int Ny, double Lx, double Ly);

/* WS_DESTRUCT - Free allocated memory and close the file (if open)
 *  ws      :   WireSlice struct
 * 
 * Returns 0 - always succeeds
 */
int ws_destruct(WireSlice_t * ws);

/* WS_READ  -   Read all data into the struct
 *  ws      :   WireSlice struct
 *  filename:   Path to the data file to open
 *  threads :   Number of subpordinate worker threads to spawn
 * 
 * Returns 0 on success
 * Returns -1 file open failure
 * Returns -2 on numerical failure
 */
int ws_read(WireSlice_t * ws, char * filename, int nthread);

/* WS_SOLVE -   Solve the problem defined by ws_init and ws_read
 *  ws      :   WireSlice struct
 *
 * Returns 0 on success
 * Returns -1 on failure
 */
int ws_solve(WireSlice_t * ws);

/* WS_WRITE -   Write the solution results to a file
 *  ws      :   WireSlice struct
 *  filename:   Path to file to write
 * 
 * Returns 0 on success
 * Returns -1 on failure
 */
int ws_write(WireSlice_t * ws, char * filename);

// The WS_READ() function uses the READ_THREAD as a worker function to
// spawn parallel threads for parsing data files.
void *read_thread(void * ws);


int ws_init(WireSlice_t * ws, unsigned int Nx, unsigned int Ny, double Lx, double Ly){
    unsigned int ii;
    
    // Dimensional parameters
    ws->Nx = Nx;    // Number of unique x wavenumbers
    ws->Ny = Ny;    // Number of unique y wavenumbers
    ws->Lx = Lx;    // Domain x-length
    ws->Ly = Ly;    // Domain y-length
    
    ws->AP = NULL;      // Upper triangular packed solution matrix
    ws->B = NULL;       // Constant vector
    ws->C = NULL;       // Coefficient vector

    ws->ndata = 0;      // Number of data points read into the problem
    ws->_halt = False;  // Used to signal threads to halt prematurely
    ws->_alock = False;  // Lock flag for multi-threading
    ws->_flock = False;
    ws->target = NULL;  // File stream for reading data
    
    // Calculate the number of coefficients
    ws->ncoef = (2*Nx+1) * (2*Ny+1);
    // Calculate the number of unique packed coefficients in AP
    // This is an upper-diagonal matrix representing an Hermitian matrix
    // so the lower diagonal can be ignored.
    ws->nAP = ws->ncoef * (ws->ncoef+1) / 2;
    
    // Allocate dynamic memory
    ws->AP = malloc(ws->nAP * sizeof(Element));
    ws->B = malloc(ws->ncoef * sizeof(Element));
    ws->C = malloc(ws->ncoef * sizeof(Element));
    if(ws->AP == NULL || ws->B == NULL || ws->C == NULL){
        ws_destruct(ws);
        return -1;
    }

    // Initialize AP and B
    for(ii=0; ii<ws->nAP; ii++)
        ws->AP[ii] = 0;
    for(ii=0; ii<ws->ncoef; ii++)
        ws->B[ii] = 0;
    // There is no need to initialize C
    
    return 0;
}


int ws_destruct(WireSlice_t * ws){
    if(ws->AP){
        free(ws->AP);
        ws->AP = NULL;
    }
    if(ws->B){
        free(ws->B);
        ws->B = NULL;
    }
    if(ws->C){
        free(ws->C);
        ws->C = NULL;
    }
    ws->ncoef=0;
    ws->ndata=0;
    ws->_halt = False;
    ws->_flock = False;
    ws->_alock = False;
    if(ws->target){
        fclose(ws->target);
        ws->target = NULL;
    }
    return 0;
}

int ws_read(WireSlice_t * ws, char *filename, int nthread){
    pthread_t threads[MAX_THREADS];
    int ii, err, thrd_err;
    // Bracket the number of threads between 1 and MAX_THREADS
    nthread = nthread <= 1 ? 1 : (nthread > MAX_THREADS ? MAX_THREADS : nthread);
    // Is a file already open!?
    if(ws->target)
        return -1;
    // Open the file
    if(!(ws->target = fopen(filename,"r")))
        return -1;
    // Initialize the halt flag and lock flags
    ws->_alock = False;
    ws->_flock = False;
    ws->_halt = False;
    // Create the specified number of threads
    err = 0;
    for(ii=0;ii<nthread;ii++)
        err = err || pthread_create(&threads[ii], NULL, &read_thread, ws);
    // If there was a failure, stop the threads and return
    if(err){
        err = -2;
        ws->_halt = True;
    }
    // Join the threads
    for(ii=0;ii<nthread;ii++){
        pthread_join(threads[ii], NULL);
    }
    // Close the file
    fclose(ws->target);
    ws->target = NULL;
    ws->_halt=False;
    return err;
}


void* read_thread(void *arg){
    WireSlice_t *ws;
    double buffer[BUFFER_SIZE * 4];
    // Intermediate matricies and vectors
    double complex  *AP,    // Upper-triangular packed Hermitian matrix
                    *B,     // Problem vector
                    *LAM;   // Coefficient vector
    // Complex floating point parameters for the wire calculation
    double complex  zt1,    // General purpose intermediate complex
                    zt2;    //
    // Floating point parameters for the wire calculation
    double          r,      // Wire radius (data file)
                    d,      // Horizontal disc location (data file)
                    theta,  // Wire angle (data file)
                    iwire,  // Total integrated wire current (data file)
                    nu_th,  // Wavenumber along the wire
                    nu_x,   // Horizontal wavenumber
                    r0,     // Radius where the wire enters the domain
                    r1,     // End radius
                    s_th,   // sine theta
                    c_th,   // cosine theta
                    ftemp;  // General purpose intermdediate double
    // Counters
    unsigned int    nread,  // Number of lines read in the current execution
                    ndata;  // Total number of lines processed in this thread
    // Indices      
    int             lam_i,  // Index in LAM
                    ap_i,   // Index in AP
                    k,      // Index in the data buffer
                    m,      // x-axis wavenumber index
                    n;      // y-axis wavenumber index
    // Helper constants
    int             izero,  // index corresponding to m=0, n=0
                    ncoef;  // coefficient relating lam_i to n

    ws = (WireSlice_t *) arg;
    
    // Initialize the intermediate matrix and vectors
    LAM = calloc(ws->ncoef, sizeof(double complex));
    B = calloc(ws->ncoef, sizeof(double complex));
    AP = calloc(ws->nAP, sizeof(double complex));


    // Calcualte the index coefficients
    // These give 
    //      lam_i = m + n*ncoef + izero
    // when m in [-Nx, Nx]
    // and  n in [-Ny, Ny]
    // and  LAM[lam_i] = Lambda(m,n)
    ncoef = (2*ws->Nx+1);
    izero = ncoef*ws->Ny + ws->Nx;

    // We'll keep track of the number of data points we've processed
    ndata = 0;
    // Continue running until one of the threads signals halt.
    while(!ws->_halt){
        // Block execution until the file lock is released
        while(ws->_flock)
            usleep(500);
        ws->_flock = True;
        nread = fread(buffer, sizeof(double), BUFFER_SIZE*4, ws->target)/4;
        ws->_flock = False;
        // If this was the last of the data, signal for the threads to halt
        if(nread < BUFFER_SIZE)
            ws->_halt = True;
        // Track the number of data points accumulated
        ndata += nread;
        // Iterate over the buffer
        for(k=0;k<4*nread;k+=4){
            r = buffer[k];
            d = buffer[k+1];
            theta = buffer[k+2];
            iwire = buffer[k+3];
            
            // calculate sine and cosine
            s_th = sin(theta);
            c_th = cos(theta);
            
            // Calculate the relevant radii
            r0 = d / c_th;
            // When calculating r1, use 1/r1 so sin(0) won't be a problem
            r1 = 1/r;
            ftemp = c_th / (d+ws->Lx);
            r1 = r1 < ftemp ? r1 : ftemp;
            ftemp = fabs(2*s_th/ws->Ly);
            r1 = r1 < ftemp ? r1 : ftemp;
            r1 = 1/r1;
            
            // Calculate Lambda
            for(n=-ws->Ny; n <= (int)ws->Ny; n++){
                for(m=-ws->Nx; m <= (int)ws->Nx; m++){
                    // Calculate the LAM index
                    lam_i = m + n*ncoef + izero;
                    // Calculate horizontal and theta wavenumbers
                    nu_x = m/ws->Lx;
                    nu_th = nu_x*c_th + n*s_th/ws->Ly;
                    // Case out zero wavenumber
                    if(nu_th == 0){
                        LAM[lam_i] = cexp(-M_2_PI*I*nu_x*d) * (r1 - r0);
                    }else{
                        zt1 = 2 * M_PI * nu_th * I;
                        zt2 = -2 * M_PI * nu_x * d * I;
                        LAM[lam_i] = (cexp(zt2 + zt1*r1) - cexp(zt2 + zt1*r0)) / zt1;
                    }
                }// n
            }// m
            // Calculate contributions to AP and B
            // repurpose m and n to be indices in Lambda
            for(m=0;m<ws->ncoef;m++){
                printf("%d: %lf + %lf j\n",m, creal(LAM[m]), cimag(LAM[n]));
                zt1 = conj(LAM[m]);
                B[m] += iwire * zt1;
                for(n=m;n<ws->ncoef;n++){
                    // The index in the packed matrix
                    ap_i = m + (n+1)*n/2;
                    AP[ap_i] += zt1 * LAM[n];
                }// n
            }// m
        }// k
    }// master while
    // Update the master struct
    // Block access until the array lock is released
    while(ws->_alock)
        usleep(500);
    ws->_alock = True;
    ws->ndata += ndata;
    ws->_alock = False;
    

    for(m=0;m<ws->ncoef;m++){
        printf("%03d: %lf + %lf j\n", m, creal(B[m]), cimag(B[m]));
    }

    printf("%08X\n%08X\n%08X\n", (int) AP, (int)B, (int)LAM);

    free(AP);
    AP = NULL;
    free(B);
    B = NULL;
    free(LAM);
    LAM = NULL;
    printf("!?\n");
    return NULL;
}
