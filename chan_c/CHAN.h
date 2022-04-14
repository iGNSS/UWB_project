#ifndef _CHAN_H
#define _CHAN_H

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

void MatrixInver(double*, int, int, double*);
void MatrixOpp(double*, int, int, double*);
double Surplus(double*, int, int);
void chan(double*, int, double*, double, double*);

#define NUM_bs 4 //使用的基站个数；

#endif
