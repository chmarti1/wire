#include <stdio.h>
#include <unistd.h>
#include "wire.h"


int main(int argc, char *argv[]){
    WireSlice_t ws;
    ws_init(&ws, 2, 2, 1, 1);
    ws_read(&ws, "test.wf", 1);
    printf("::%d\n", ws.ndata);
    ws_destruct(&ws);
    return 0;
}



