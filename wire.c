#include <stdio.h>
#include <unistd.h>
#include "wire.h"


int main(int argc, char *argv[]){
    WireSlice_t ws;
    int ii, err;
    err = ws_init(&ws, 10, 10, 1, 1);
    if(!err){
        printf("Reading file...\n");
        ws_read(&ws, "test.wf", 1);
    }
    printf("::%d\n", ws.ndata);
    ws_solve(&ws);
    ws_write(&ws, "output.wf");
    ws_destruct(&ws);
    return 0;
}



