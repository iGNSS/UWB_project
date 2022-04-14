#include "CHAN.h"

int main() {
    double *bs = NULL;
    double *dis = NULL;
    double *dis1 = NULL;
    double *dis2 = NULL;
    double *theta = NULL;
    int bs_row = NUM_bs;    //使用的基站个数，需作为输入传入chan函数；
    double sigma = 10;       //基站测量误差范围；

    bs = (double*)malloc(sizeof(double)*NUM_bs*2);  //基站坐标数组{x0,y0,x1,y1,...,x5,y5};
    dis = (double*)malloc(sizeof(double)*NUM_bs);   //基站距离数组{d0,d1,...,d5};
    dis1 = (double*)malloc(sizeof(double)*NUM_bs);   //基站距离数组{d0,d1,...,d5};
    dis2 = (double*)malloc(sizeof(double)*NUM_bs);   //基站距离数组{d0,d1,...,d5};
    theta = (double*)malloc(sizeof(double)*2);      //数字钥匙坐标{x,y};

    bs[0] = -161;             //A 
    bs[1] = -83;              
    bs[2] = 253;              //B
    bs[3] = -80.5;
    bs[4] = 0;                //E
    bs[5] = 0;
    bs[6] = 117;              //F
    bs[7] = 0;
    bs[8] = 249;              //C
    bs[9] = 84; 
    bs[10] = -159;            //D
    bs[11] = 89;

    dis[0] = 1263;
    dis[1] = 910;
    dis[2] = 1139;
    dis[3] = 1056;
    dis[4] = 1378;
    dis[5] = 1020;
    
    
    dis1[0] = 1019;
    dis1[1] = 626;
    dis1[2] = 892;
    dis1[3] = 819;
    dis1[4] = 1150;
    dis1[5] = 729;
    
    dis2[0] = 478;
    dis2[1] = 347;
    dis2[2] = 415;
    dis2[3] = 402;
    dis2[4] = 538;
    dis2[5] = 609;

    chan(bs, bs_row, dis, sigma, theta);  //调用chan函数；
    printf("%f,%f\n",theta[0],theta[1]);
    chan(bs, bs_row, dis1, sigma, theta);  //调用chan函数；
    printf("%f,%f\n",theta[0],theta[1]);
    chan(bs, bs_row, dis2, sigma, theta);  //调用chan函数；
    printf("%f,%f\n",theta[0],theta[1]);

    getchar();
    return 0;
}