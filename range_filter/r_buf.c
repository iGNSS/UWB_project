#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#define BUF_DEPTH 7

buf0 = (double*)malloc(sizeof(double)*BUF_DEPTH);
buf1 = (double*)malloc(sizeof(double)*BUF_DEPTH);
buf2 = (double*)malloc(sizeof(double)*BUF_DEPTH);
buf3 = (double*)malloc(sizeof(double)*BUF_DEPTH);
buf4 = (double*)malloc(sizeof(double)*BUF_DEPTH);
buf5 = (double*)malloc(sizeof(double)*BUF_DEPTH);

void r_buf(double arr, double* buf0, double* buf1, double* buf2, double* buf3, double* buf4, double* buf5){
    int i;
      
    for(i=0,i<BUF_DEPTH-1,i++){
        buf0[i+1] = buf0[i];
    }
    buf0[0] = arr[0];

    for(i=0,i<BUF_DEPTH-1,i++){
        buf1[i+1] = buf1[i];
    }
    buf1[0] = arr[1];

    for(i=0,i<BUF_DEPTH-1,i++){
        buf2[i+1] = buf2[i];
    }
    buf2[0] = arr[2];

    for(i=0,i<BUF_DEPTH-1,i++){
        buf3[i+1] = buf3[i];
    }
    buf3[0] = arr[3];

    for(i=0,i<BUF_DEPTH-1,i++){
        buf4[i+1] = buf4[i];
    }
    buf4[0] = arr[4];

    for(i=0,i<BUF_DEPTH-1,i++){
        buf5[i+1] = buf5[i];
    }
    buf5[0] = arr[5];

}