#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#define ERROR 0
#define DEVIATION 2

//range_arr_d = range_arr;

void range_filter(double* range_mat, double* range_mat_d, int arr_col, int buf_size, double* range_buf){
    int i,j;
    double range_arr[6];
    double range_buf_d[7][6];

    for (i=0;i<arr_col;i++){
        range_arr[i] = *(range_mat+i);
    }

    //'0' can not appear in first range mat;
    for (i=0;i<arr_col;i++){
        if (*(range_arr_d+i)==0){
            return ERROR;
        }
    }

    //fix '0' in range arr;
    for (i=0;i<arr_col;i++){
        if (range_arr[i]==0){
            range_arr[i] = *(range_arr_d+i);
        }
    }

    //fix BS deviation;
    for (i=0;i<arr_col;i++){
        range_arr[i] = range_arr[i]-DEVIATION;
        if (range_arr[i]<=0){
            range_arr[i] = 0;
        }
    }

    //rang buffer reg;
    for (i=0;i<buf_size;i++){
        for (j=0;j<arr_col;j++){
            range_buf_d[i][j] = *(range_buf+i*arr_col+j);
        }
    }

    //buffer first in first out;
    for (j=0;j<arr_col;j++){
        for (i=0;i<buf_size;i++){
            range_buf_d[i+1][j] = range_buf_d[i][j];
        }
    }
    for (j=0;j<arr_col;j++){
        range_buf_d[0][j] = range_arr[j];
    }

    
    
    








    
}


/*----------------------------------------------------------------------------------------*/
/*-------------------------------------SG_FILTER------------------------------------------*/
/*----------------------------------------------------------------------------------------*/

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
