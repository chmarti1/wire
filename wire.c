#include <stdio.h>
#include <unistd.h>
#include "wire.h"


int main(int argc, char *argv[]){
    WireSlice_t ws;
    int m,n,ii, err;
    err = ws_init(&ws, 1, 1, 1, 1);
    if(!err){
        printf("Reading file...\n");
        ws_read(&ws, "test.wf", 8);
    }
    printf("AP = [\n");
    for(m=0;m<ws.ncoef;m++){
        printf("[");
        for(n=0;n<ws.ncoef;n++){
            if(m<=n){
                ii = m + n*(n+1)/2;
                printf("%6.1lf+%6.1lfj,  ", creal(ws.AP[ii]), cimag(ws.AP[ii]));
            }else{
                ii = n + m*(m+1)/2;
                printf("%6.1lf+%6.1lfj,  ", creal(ws.AP[ii]), cimag(ws.AP[ii]));
            }
        }
        printf("],\n");
    }
    printf("]\nB = [\n");
    for(m=0;m<ws.ncoef;m++){
        printf("%6.1lf+%6.1lfj,\n", creal(ws.B[m]), cimag(ws.B[m]));
    }
    printf("]\nC = [\n");
    ws_solve(&ws);
    ws_write(&ws, "output.wf");
    
    for(m=0;m<ws.ncoef;m++){
        printf("%7.4lf+%7.4lfj,\n", creal(ws.B[m]), cimag(ws.B[m]));
    }
    printf("]\n");
    ws_destruct(&ws);
    return 0;
}



