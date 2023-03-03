#include <complex.h>
#include <cblas.h>
#include <lapacke.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>
#include <time.h>

// XXX - may need the -fopenmp compiler directive for parallel processing? 
// XXX - Need "-lpthread -lm -llapacke" in the linking stage

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
    unsigned int Nx, Ny;
    double Lx, Ly;
    double dshift;
    unsigned int ncoef;
    unsigned int nAP;
    unsigned int ntotal, ndata, nused;
    unsigned char _halt;
    unsigned char _verbose;
    time_t _start;
    pthread_mutex_t _matlock;
    pthread_mutex_t _filelock;
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
 *  verbose :   Should the operations print status messages?
 * 
 * Returns 0 on success
 * Returns -1 on memory allocation failure
 */
int ws_init(WireSlice_t * ws, unsigned int Nx, unsigned int Ny, double Lx, double Ly, unsigned char verbose);

/* WS_DESTRUCT - Free allocated memory and close the file (if open)
 *  ws      :   WireSlice struct
 * 
 * Returns 0 - always succeeds
 */
int ws_destruct(WireSlice_t * ws);

/* WS_VERBOSE - Turn stdout status messages on or off
 *  ws      :   WireSlice struct
 *  verbose :   0 or 1
 * 
 * Returns 0 always
 */
 int ws_verbose(WireSlice_t * ws, unsigned char verbose);

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

// The EJTH() function is about twice as fast as calling cexp() for
// evaluating exp(j*theta).  Because there is no real component, it
// can work with a single call to cosine instead of a power call.
double complex ejtheta(double theta);

int ws_init(WireSlice_t * ws, unsigned int Nx, unsigned int Ny, double Lx, double Ly, unsigned char verbose){
    int err = 0;
    
    // Set the time of the algorithm start
    ws->_start = time(NULL);
    
    // Dimensional parameters
    ws->Nx = Nx;    // Number of unique x wavenumbers
    ws->Ny = Ny;    // Number of unique y wavenumbers
    ws->Lx = Lx;    // Domain x-length
    ws->Ly = Ly;    // Domain y-length
    ws->dshift = 0.;
    
    ws->AP = NULL;      // Upper triangular packed solution matrix
    ws->B = NULL;       // Solution vector

    ws->ntotal = 0;     // Total number of data point available to ws_read()
    ws->ndata = 0;      // Number of data points found by ws_read()
    ws->nused = 0;      // Number of data points that actually intersect the domain
    ws->_halt = False;  // Used to signal threads to halt prematurely
    ws->_verbose = verbose;
    // Initialize the mutexes for locking the master matrix and file
    // Default shared state for pthread mutexes is PRIVATE - only one process.
    if(     pthread_mutex_init(&ws->_matlock, NULL) ||
            pthread_mutex_init(&ws->_filelock, NULL) ){
        fprintf(stderr, "WS_INIT: Failed while initializing mutexes!?\n");
        ws_destruct(ws);
        return -1;
    }
    
    ws->target = NULL;  // File stream for reading data
    
    // Calculate the number of coefficients
    ws->ncoef = (2*Nx+1) * (2*Ny+1) + 1;
    // Calculate the number of unique packed coefficients in AP
    // This is an upper-diagonal matrix representing an Hermitian matrix
    // so the lower diagonal can be ignored.
    ws->nAP = ws->ncoef * (ws->ncoef+1) / 2;
    
    // Allocate dynamic memory
    ws->AP = calloc(ws->nAP, sizeof(Element));
    ws->B = calloc(ws->ncoef, sizeof(Element));
    if(ws->AP == NULL || ws->B == NULL){
        fprintf(stderr, "WS_INIT: Failed during memory allocation, nAP=%d, ncoef=%d\n",ws->nAP, ws->ncoef);
        ws_destruct(ws);
        err = -1;
    }

    if(ws->_verbose)
        printf("[%6.0f] WS_INIT RC: %d\n", difftime(time(NULL), ws->_start), err);    
    return err;
}





int ws_destruct(WireSlice_t * ws){
    int err = 0;
    
    if(ws->AP){
        free(ws->AP);
        ws->AP = NULL;
    }
    if(ws->B){
        free(ws->B);
        ws->B = NULL;
    }

    ws->ncoef=0;
    ws->ndata=0;
    ws->_halt = False;

    pthread_mutex_destroy(&ws->_matlock);
    pthread_mutex_destroy(&ws->_filelock);

    if(ws->target){
        fclose(ws->target);
        ws->target = NULL;
    }
    
    if(ws->_verbose)
        printf("[%6.0f] WS_DESTRUCT RC: %d\n", difftime(time(NULL), ws->_start), err);
    return err;
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
    // Detect the nuber of available data points
    fseek(ws->target, 0, SEEK_END);
    ws->ntotal = ftell(ws->target) / (4 * sizeof(double));
    fseek(ws->target, 0, SEEK_SET);
    if(ws->_verbose){
        printf("  Using %d threads to read from file: %s\n", nthread, filename);
        printf("    Nx=%d, Ny=%d\n", ws->Nx, ws->Ny);
        printf("    Lx=%lf, Ly=%lf\n", ws->Lx, ws->Ly);
    }
    
    // Set the threads to run until halted by one of them
    ws->_halt = False;
    // Create the specified number of threads
    err = 0;
    for(ii=0;ii<nthread;ii++)
        err = err || pthread_create(&threads[ii], NULL, &read_thread, ws);
    // If there was a failure, stop the threads and return
    if(err){
        err = -2;
        ws->_halt = True;
        fprintf(stderr,"WS_READ: encountered an error during thread creation.\n");
    }
    // Join the threads
    for(ii=0;ii<nthread;ii++){
        pthread_join(threads[ii], NULL);
    }
        
    // Close the file
    fclose(ws->target);
    ws->target = NULL;
    ws->_halt=False;

    // Add a newline to the front of the RC status to close out the
    // threads' status messages
    if(ws->_verbose)
        printf("\n[%6.0f] WS_READ RC: %d\n", difftime(time(NULL), ws->_start), err);
    return err;
}


int ws_solve(WireSlice_t *ws){
    lapack_int err;
    //lapack_int *ipiv;
    //ipiv = malloc(ws->ncoef * sizeof(lapack_int));
    
    /* Prototype from LAPACKE for Hermitian matrix
    lapack_int LAPACKE_zppsv( int matrix_layout, char uplo, lapack_int n,
                          lapack_int nrhs, lapack_complex_double* ap,
                          lapack_complex_double* b, lapack_int ldb );
    */
    
    err = LAPACKE_zppsv(    LAPACK_COL_MAJOR, 'U', ws->ncoef, 1, 
                            ws->AP, ws->B, ws->ncoef);
    
    
    
    /* Prototype from LAPACKE for symmetrical complex matrix
    lapack_int LAPACKE_zspsv( int matrix_layout, char uplo, lapack_int n,
                          lapack_int nrhs, lapack_complex_double* ap,
                          lapack_int* ipiv, lapack_complex_double* b,
                          lapack_int ldb );
    */
    /*
    err = LAPACKE_zspsv(    LAPACK_COL_MAJOR, 'U', ws->ncoef, 1, 
                            ws->AP, ipiv, ws->B, ws->ncoef);
    free(ipiv);
    */
    
    if(ws->_verbose)
        printf("[%6.0f] WS_SOLVE RC: %d\n", difftime(time(NULL), ws->_start), err);
    return err;
}


int ws_write(WireSlice_t *ws, char * filename){
    FILE * fd = NULL;
    int err = 0;
    
    fd = fopen(filename,"w");
    
    if(fd){
        // Write header information
        fwrite(&ws->Nx, sizeof(unsigned int), 1, fd);
        fwrite(&ws->Ny, sizeof(unsigned int), 1, fd);
        fwrite(&ws->Lx, sizeof(double), 1, fd);
        fwrite(&ws->Ly, sizeof(double), 1, fd);
        fwrite(&ws->dshift, sizeof(double), 1, fd);
        
        // Write complex vector solution
        fwrite(ws->B, sizeof(double complex), ws->ncoef, fd);
        
        fclose(fd);
        fd = NULL;
    }else{
        err = -1;
        fprintf(stderr,"WS_WRITE: Failed to open file for writing:\n");
        fprintf(stderr,"  %s", filename);
    }
    
    if(ws->_verbose)
        printf("[%6.0f] WS_WRITE RC: %d\n", difftime(time(NULL), ws->_start), err);
    return err;
}

/*
 *  HELPER FUNCTION DEFINITIONS
 */
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
                    h0,     // height for domain intersection test
                    ft1,    // General purpose intermdediate double
                    ft2;    //
    // Counters
    unsigned int    ndata_this,  // Number of lines read in the current execution
                    nused_this,  // Number of points used in the current exection
                    ndata,  // Total number of lines read by this thread
                    nused;  // Total number of points used by this thread

    // Indices      
    int             lam_i,  // Index in LAM
                    ap_i,   // Index in AP
                    k,      // Index in the data buffer
                    m,      // x-axis wavenumber index
                    n;      // y-axis wavenumber index
    // Helper constants
    int             izero,  // index corresponding to m=0, n=0
                    icoef;  // coefficient relating lam_i to n

    ws = (WireSlice_t *) arg;
    
    // Initialize the intermediate matrix and vectors
    LAM = calloc(ws->ncoef, sizeof(double complex));
    B = calloc(ws->ncoef, sizeof(double complex));
    AP = calloc(ws->nAP, sizeof(double complex));

    // Calcualte the index coefficients
    // These give the location in the Lambda vector for given vertical,
    // horizontal wavenumber pair, (m,n)
    //      lam_i = m + n*icoef + izero
    // when m in [-Nx, Nx]
    // and  n in [-Ny, Ny]
    // and  LAM[lam_i] = Lambda(m,n)
    icoef = (2*ws->Nx+1);
    izero = icoef*ws->Ny + ws->Nx;

    // We'll keep track of the number of data points we've processed
    ndata = 0;
    nused = 0;
    // Continue running until one of the threads signals halt.
    while(!ws->_halt){
        // Block execution until the file lock is released
        pthread_mutex_lock(&ws->_filelock);
        ndata_this = fread(buffer, sizeof(double), BUFFER_SIZE*4, ws->target)/4;
        pthread_mutex_unlock(&ws->_filelock);
        // If this was the last of the data, signal for the threads to halt
        if(ndata_this < BUFFER_SIZE)
            ws->_halt = True;
        // Track the number of data points accumulated
        ndata += ndata_this;
        nused_this = 0;
        // Iterate over the buffer
        for(k=0;k<4*ndata_this;k+=4){
            r = buffer[k];
            d = buffer[k+1] + ws->dshift;
            theta = buffer[k+2];
            iwire = buffer[k+3];
            
            // calculate sine and cosine
            zt1 = ejtheta(theta);
            c_th = creal(zt1);
            s_th = cimag(zt1);
            
            // Calculate the relevant radii
            // r0 is the projection along the wire to the first domain edge
            r0 = d / c_th;
            // Calculate the height at which the wire intersects the domain
            // This isn't actually important - it's just the easiest way
            // to verify that the wire passes through the domain.
            h0 = fabs(r0 * s_th);
            // r1 is the smallest of
            //  (1) the wire radius
            //  (2) the wire projection to the far edge of the domain
            //  (3) the wire projection to either of the horizontal 
            //      domain boundaries (top or bottom). 
            // When the wire is horizontal, (3) is infinite, so we first
            // calculate 1/r1 and accept the largest value.
            // When the wire terminates inside the domain r1 = r
            // When the wire doesn't reach the domain, r1 < r0
            r1 = 1/r;
            ft1 = c_th / (d+ws->Lx);    // Candidate: r1=(d+Lx)/cos(theta)
            r1 = r1 > ft1 ? r1 : ft1;
            ft1 = fabs(2*s_th/ws->Ly);  // Candidate: r1=Ly/2sin(theta)
            r1 = r1 > ft1 ? r1 : ft1;
            r1 = 1/r1;
            // Only proceed if the wire intersects the domain
            // There are two ways for the wire to miss:
            //  (1) the angle is too wide
            //  (2) the disc is too far away and the wire doesn't reach
            if(h0 < ws->Ly/2 && r1 > r0){
                nused_this ++;
                // Calculate the Lambda vector
                // The reverse of Lambda is its complex conjugate, so
                // only half of the coefficients need to be calculated.
                // We'll iterate over ALL n-values, but only half of the
                // m-values.
                for(n=-ws->Ny; n <= (int)ws->Ny; n++){
                    // Deal with m=0 specially
                    lam_i = n*icoef + izero;
                    //nu_x = 0;     No need to calculate nu_x
                    nu_th = n*s_th/ws->Ly;
                    // Case out zero wavenumber
                    if(n == 0 || nu_th ==0){
                        LAM[lam_i] = (r1 - r0);
                    }else{
                        ft1 = 2 * M_PI * nu_th;
                        LAM[lam_i] = (ejtheta(ft1*r1) - ejtheta(ft1*r0)) / (ft1 * I);
                    }
                    // Deal with m!=0 in complex conjugate pairs
                    for(m=1; m <= (int)ws->Nx; m++){
                        // Calculate the LAM index
                        lam_i = m + n*icoef + izero;
                        // Calculate horizontal and theta wavenumbers
                        nu_x = m/ws->Lx;
                        nu_th = nu_x*c_th + n*s_th/ws->Ly;
                        // Phase shift implied by the disc location
                        ft2 = -2 * M_PI * nu_x * d;
                        // Case out zero wavenumber
                        if(nu_th == 0){
                            zt1 = ejtheta(ft2) * (r1-r0);
                        }else{
                            ft1 = 2 * M_PI * nu_th;
                            zt1 = (ejtheta(ft1*r1 + ft2) - ejtheta(ft1*r0 + ft2)) / (ft1 * I);
                        }
                        LAM[lam_i] = zt1;
                        // Assign the complex conjugate
                        lam_i = -m -n*icoef + izero;
                        LAM[lam_i] = conj(zt1);
                    }// m
                }// n
                // Assign the offset current coefficient
                LAM[ws->ncoef-1] = 1.;
                
                // Calculate contributions to AP and B
                // repurpose m and n to be indices in Lambda
                for(m=0;m<ws->ncoef;m++){
                    // Take the conjugate 
                    // AP = sum LAM* LAM
                    //  B = sum LAM* Iwire
                    zt1 = conj(LAM[m]);
                    //zt1 = LAM[m];
                    B[m] += iwire * zt1;
                    for(n=m;n<ws->ncoef;n++){
                        // The index in the packed matrix
                        ap_i = m + (n+1)*n/2;
                        AP[ap_i] += zt1 * LAM[n];
                    }// n
                }// m
            }// intersection test
        }// k
        // Update the data counters
        nused += nused_this;
        pthread_mutex_lock(&ws->_matlock);
        ws->ndata += ndata_this;
        ws->nused += nused_this;
        pthread_mutex_unlock(&ws->_matlock);
        // Update status information
        if(ws->_verbose){
            fprintf(stdout, "\x1B""[1K""\x1B""[1G""Status: %d Total / %d(%4.1f%%) Processed / %d(%4.1f%%) Used", 
                ws->ntotal, 
                ws->ndata, 100*((float) ws->ndata)/ws->ntotal, 
                ws->nused, 100*((float) ws->nused)/ws->ndata);
            fflush(stdout);
        }
    }// master while
    // Update the master struct
    // Block access until the matrix lock is released
    pthread_mutex_lock(&ws->_matlock);
    for(m=0;m<ws->ncoef;m++){
        ws->B[m] += B[m];
    }
    for(m=0;m<ws->nAP;m++){
        ws->AP[m] += AP[m];
    }
    pthread_mutex_unlock(&ws->_matlock);


    free(AP);
    AP = NULL;
    free(B);
    B = NULL;
    free(LAM);
    LAM = NULL;
    return NULL;
}

// calculate exp(j*theta) efficiently
// This is approximately 42ms per execution on jane, while
// cexp( ) is approimately 66ms per execution.
double complex ejtheta(double theta){
    double c_th, mod_theta;
    c_th = cos(theta);
    mod_theta = theta / (2*M_PI);
    mod_theta -= floor(mod_theta);
    if(mod_theta >= 0.5)
        return CMPLX(c_th, -sqrt(1. - c_th*c_th));
    return CMPLX(c_th, sqrt(1. - c_th*c_th));
}
