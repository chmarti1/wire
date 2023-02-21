#include <stdio.h>
#include <unistd.h>
#include "wire.h"


int main(int argc, char *argv[]){
    WireSlice_t ws;
    int m,n,ii, err;
    err = ws_init(&ws, 1, 1, 1, 1);
    if(!err){
        printf("Reading file...\n");
        ws_read(&ws, "test.wf", 1);
    }
    printf("::%d\n", ws.ndata);
    printf("AP\n");
    for(m=0;m<ws.ncoef;m++){
        printf("( ");
        for(n=0;n<ws.ncoef;n++){
            if(m<=n){
                ii = m + n*(n+1)/2;
                printf("%6.1lf+%6.1lfj  ", creal(ws.AP[ii]), cimag(ws.AP[ii]));
            }else{
                printf("      *         ");
            }
        }
        printf(")  ( %6.1lf+%6.1lfj )\n", creal(ws.B[m]), cimag(ws.B[m]));
    }
    ws_solve(&ws);
    ws_write(&ws, "output.wf");
    ws_destruct(&ws);
    return 0;
}



