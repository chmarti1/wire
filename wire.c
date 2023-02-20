#include <stdio.h>
#include <unistd.h>
#include "wire.h"


int main(int argc, char *argv[]){
    WireSlice_t ws;
    ws_init(&ws, 10, 10, 1, 1);
    ws_read(&ws, "test.wf", 8);
    printf("::%d\n", ws.ndata);
    ws_solve(&ws);
    ws_destruct(&ws);
    return 0;
}



