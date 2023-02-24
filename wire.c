#include <stdio.h>
#include <unistd.h>
#include "wireslice.h"


const char [] help_text = "wire -h <filename>

int main(int argc, char *argv[]){
    WireSlice_t ws;
    int err;
    err = ws_init(&ws, 20, 20, 1, 1);
    err = err || ws_read(&ws, "test.wf", 8);
    err = err || ws_solve(&ws);
    err = err || ws_write(&ws, "output.wf");
    ws_destruct(&ws);
    return err;
}



