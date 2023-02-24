#include <stdio.h>
#include <unistd.h>
#include "wire.h"


int main(int argc, char *argv[]){
    WireSlice_t ws;
    int m,n,ii, err;
    err = ws_init(&ws, 20, 20, 1, 1);
    if(!err){
        printf("Reading file...\n");
        ws_read(&ws, "test.wf", 4);
    }
    printf("Solving...\n");
    ws_solve(&ws);
    ws_write(&ws, "output.wf");
    
    ws_destruct(&ws);
    return 0;
}



