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
#define TWOPI       2 * M_PI    // Duh
#define NREAD       5           // How many floating point numbers to read?


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
 * The binary file should contain groups of five double-precision values
 * in order: (R, X, Y, THETA, I).  R is the wire radius, X,Y is the 
 * location of the disc center, THETA is the wire angle, and I is the 
 * measured wire current. 
 * 
 * WS_READ() calls the READ_THREAD() function in a number of parallel
 * threads to stream in data from the source file into the WireSlice
 * struct's A matrix and B vector.  
 * 
 * After WS_SOLVE() returns successfully, the solution vector can be 
 * written to a file using WS_WRITE().
 * 
 */

typedef struct WireSlice {
    Element * AP;
    Element * B;
    unsigned int Nx, Ny;
    double Lx, Ly;
    double xshift, yshift;
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
 * Returns the return value of the LAPACK's ZPPSV function. It is 
 * negative if there is an error in the configuration parameters, and 
 * it is positive if the matrix is singular or there is some other error
 * in the process.  The return value is zero on success.
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

/* WS_SHIFT -   Set the x- and y- shift values to apply to all data
 *  ws      :   WireSlice struct
 *  x,y     :   x- and y-shift values
 * 
 * Returns 0
 */
int ws_shift(WireSlice_t * ws, double x, double y);

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
    ws->xshift = 0.;
    ws->yshift = 0.;
    
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
    printf("nthread: %d\n", nthread);
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


int ws_shift(WireSlice_t * ws, double x, double y){
    ws->xshift = x;
    ws->yshift = y;
    if(ws->_verbose)
        printf("[%6.0f] WS_SHIFT RC: 0   (%lf, %lf)\n", difftime(time(NULL), ws->_start), x, y);
    return 0;
}

/*
 *  HELPER FUNCTION DEFINITIONS
 */

void* read_thread(void *arg){
    WireSlice_t *ws;
    double buffer[BUFFER_SIZE * NREAD];
    // Intermediate matricies and vectors
    double complex  *AP,    // Upper-triangular packed Hermitian matrix
                    *B,     // Problem vector
                    *LAM;   // Coefficient vector
    // Complex floating point parameters for the wire calculation
    double complex  zt1,    // General purpose intermediate complex
                    zt2;    //
    // Floating point parameters for the wire calculation
    double          r,      // Wire radius (data file)
                    r0,     // Radius where the wire enters the domain
                    r1,     // End radius
                    x0,     // Horizontal disc location (data file)
                    y0,     // Vertical disc location (data file)
                    theta,  // Wire angle (data file)
                    iwire,  // Total integrated wire current (data file)
                    nu_th,  // Wavenumber along the wire
                    phi,  // Phase induced by the disc center location
                    nu_x,   // X- and Y- components of the wavenumber
                    nu_y,   //   vector.
                    s_th,   // sine theta
                    c_th,   // cosine theta
                    ft1,    // General purpose intermdediate double
                    ft2;    //
    // Registers for calculating intersection
    double          rho0, 
                    rho1,
                    rho,
                    d;
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
        ndata_this = fread(buffer, sizeof(double), BUFFER_SIZE*5, ws->target)/NREAD;
        pthread_mutex_unlock(&ws->_filelock);
        // If this was the last of the data, signal for the threads to halt
        if(ndata_this < BUFFER_SIZE)
            ws->_halt = True;
        // Track the number of data points accumulated
        ndata += ndata_this;
        nused_this = 0;
        // Iterate over the buffer
        for(k=0;k<NREAD*ndata_this;k+=NREAD){
            // Extract the five values from each entry
            r = buffer[k];
            x0 = buffer[k+1] + ws->xshift;
            y0 = buffer[k+2] + ws->yshift;
            theta = buffer[k+3];
            iwire = buffer[k+4];
            
            // calculate sine and cosine of the wire trajectory
            zt1 = ejtheta(theta);
            c_th = creal(zt1);
            s_th = cimag(zt1);
            
            /* Calculate wire intersection with the bounding edges. We
             * will name the edges AB, BC, CD, and DA
             *         ^
             *         |            +  (x0,y0)
             *   B ,---|---, C
             *     |   |   |
             * --------+------------>
             *     |   |   |
             *   A '---|---' D
             * 
             * A line projected from x0,y0 with radius r and angle theta
             * may pass through 0, 1, or 2 of the domain edges.  This 
             * algorithm detects the maximum and minimum radii of the
             * wire's path through the domain (if any).
             * 
             * The algorithm works in rho=1/r instead of r to avoid 
             * singularities when the line is parallel with an edge.
             * There is a pair of sorted values, rho0 and rho1, that
             * indicate the maximum and minimum values of rho at an 
             * edge intersection.  Initially, rho0 == rho1 == 1/r
             * so that if the intersection radius is beyond the wire tip,
             * rho < rho0, and the possible intersection will be ignored.
             * This has the double benefit of ignoring singular cases
             * and edges outside of the wire's reach in a single 
             * comparison.
             * 
             * Then, the algorithm calculates the intersection location, 
             * d, along the edge from the nearest axis.  If it lies 
             * within the length of the edge segment, the wire 
             * intersects the edge, and the value of rho is inserted
             * into the rho0, rho1 pair.
             * 
             * After all four edges have been tested for intersection,
             * the rho0 and rho1 pair are tested for equality.  If
             * rho1 > rho0, there was an intersection.  If only one edge
             * was intersected, rho0 will still be 1/r (the wire tip).
             * If two intersections were found, then rho0 and rho1 will
             * be appropriate radii for the entrance and exit of the 
             * wire.
             */
             
            rho0 = 1/r;
            rho1 = rho0;
            
            // Start with AB
            rho = -c_th / (x0 + ws->Lx/2);
            // Only bother if rho is greater than the smallest candidate
            if(rho > rho0){
                d = s_th / rho + y0;
                if(fabs(d) <= ws->Ly/2){
                    if(rho>rho1){
                        rho0 = rho1;
                        rho1 = rho;
                    }else{
                        rho0 = rho;
                    }
                }
            }
            
            // next, check CD
            rho = c_th / (-x0 + ws->Lx/2);
            // Only bother if rho is greater than the smallest candidate
            if(rho > rho0){
                d = s_th / rho + y0;
                if(fabs(d) <= ws->Ly/2){
                    if(rho>rho1){
                        rho0 = rho1;
                        rho1 = rho;
                    }else{
                        rho0 = rho;
                    }
                }
            }
            
            // next, check BC
            rho = s_th / (-y0 + ws->Ly/2);
            // Only bother if rho is greater than the smallest candidate
            if(rho > rho0){
                d = c_th / rho + x0;
                if(fabs(d) < ws->Lx/2){
                    if(rho>rho1){
                        rho0 = rho1;
                        rho1 = rho;
                    }else{
                        rho0 = rho;
                    }
                }
            }
            
            // next, check DA
            rho = -s_th / (y0 + ws->Ly/2);
            // Only bother if the candidate is greater than the smallest candidate
            if(rho > rho0){
                d = c_th / rho + x0;
                if(fabs(d) < ws->Lx/2){
                    if(rho>rho1){
                        rho0 = rho1;
                        rho1 = rho;
                    }else{
                        rho0 = rho;
                    }
                }
            }
            
            // If the wire intersects the domain
            if(rho1>rho0){
                // Convert rhos back into radii
                r1 = 1./rho0;
                r0 = 1./rho1;
                
                nused_this ++;
                // Calculate the Lambda vector
                // The reverse of Lambda is its complex conjugate, so
                // only half of the coefficients need to be calculated.
                // We'll iterate over ALL n-values, but only half of the
                // m-values.
                for(n=-ws->Ny; n <= (int)ws->Ny; n++){
                    // Calculate the vertical component of wavenumber
                    nu_y = TWOPI * n/ws->Ly;
                    //
                    // Deal with m=0 specially, outside of the loop
                    // below.
                    lam_i = n*icoef + izero;
                    // calculate the phase and wavenumber along the wire
                    // When m==0, only the vertical component matters
                    phi = nu_y * y0;
                    nu_th = nu_y * s_th;
                    // Catch the special case of a zero wavenumber
                    if(n == 0 || nu_th ==0){
                        LAM[lam_i] = (r1 - r0) * ejtheta(phi);
                    }else{ 
                        LAM[lam_i] = (ejtheta(nu_th*r1 + phi) - ejtheta(nu_th*r0 + phi)) / (nu_th * I);
                    }
                    //
                    // Now, deal with m!=0 in complex conjugate pairs
                    // in a loop.  We'll index over ONLY the positive
                    // values of m, then negate BOTH m and n.
                    for(m=1; m <= (int)ws->Nx; m++){
                        // Calculate the horizontal wavenumber component
                        nu_x = TWOPI * m/ws->Lx;
                        // Calculate the LAM index
                        lam_i = m + n*icoef + izero;
                        // This is actually 2*pi*phi
                        phi = nu_x*x0 + nu_y*y0;
                        nu_th = nu_x*c_th + nu_y*s_th;
                        // Case out zero wavenumber
                        if(nu_th == 0){
                            zt1 = (r1-r0) * ejtheta(phi);
                        }else{
                            zt1 = (ejtheta(nu_th*r1 + phi) - ejtheta(nu_th*r0 + phi)) / (nu_th * I);
                        }
                        LAM[lam_i] = zt1;
                        // Assign the complex conjugate at -m,-n
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

