#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#define N 0             //System precision;
#define BUF_DEPTH 7     //Buffer depth, number must be odd;

double  sgfilter(double*, int);

//function instance; 
int main(){
    double arr[BUF_DEPTH] = {2,1,2.2,4,4.1,9.6,11};
    int size_arr = BUF_DEPTH;
    double dout;

    dout = sgfilter(arr, size_arr);

    printf("%f\n",dout);
    return 0;
}

double sgfilter(double arr[], int size_arr){
    double index;
    int i;
    double range_sg;                            //filtered value, delay is (BUF_DEPTH-1)/2;
    double x[size_arr];
    double a,b,c,m1,m2,m3,z1,z2,z3;a=b=c=0;
    double sumx=0,sumx2=0,sumx3=0,sumx4=0,sumy=0,sumxy=0,sumx2y=0;

    for(i=0;i<size_arr;i++){
        x[i] = i;
    }

    for(i=0;i<size_arr;i++){
        sumx+=x[i];sumy+=arr[i];
        sumx2+=pow (x[i],2); sumxy+=x[i]*arr[i];
        sumx3+=pow(x[i],3); sumx2y+=pow(x[i],2)*arr[i];
        sumx4+=pow(x[i],4);
    }

    do{
        m1=a; 
        a=(sumx2y-sumx3*b-sumx2*c)/sumx4; 
        z1=(a-m1)*(a-m1);
        m2=b; 
        b=(sumxy-sumx*c-sumx3*a)/sumx2; 
        z2=(b-m2)*(b-m2);
        m3=c; 
        c=(sumy-sumx2*a-sumx*b)/size_arr; 
        z3=(c-m3)*(c-m3);
    }while((z1>N)||(z2>N)||(z3>N));

    index = (size_arr/2);   
    range_sg = (a*index*index + b*index + c);
    //printf ("y=%9.6fx*x+%9.6fx+%9.6f\n",a,b,c);
    return range_sg;
}