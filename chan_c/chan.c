#include "CHAN.h"

/*----------------------------------------------------------------------------------------*/
/*---------------------------------------CHAN---------------------------------------------*/
/*----------------------------------------------------------------------------------------*/
void chan(double* bs,int bs_row,double* dis,double sigma,double* theta){
    int i = 0;
    int j = 0;
    int k = 0;
    int row_col = 0;
    double sum = 0; 
    double *arr = NULL;
    double *result = NULL;
    double bs_mat[bs_row][2];
    double k1[bs_row];
    double g1[bs_row][3];
    double g1_t[3][bs_row];
    double h1[bs_row];
    double mat_temp0[3][bs_row];
    double mat_temp1[3][3];
    double mat_temp1_inv[3][3];
    double theta0[2] = {0,0};
    double dis_t[bs_row];
    double b1[bs_row][bs_row];  
    double cov1[bs_row][bs_row];
    double cov1_inv[bs_row][bs_row];
    double mat_temp2[3][bs_row];
    double theta1[3] = {0,0,0};
    double h2[3] = {0,0,0};
    double g2[3][2] = {{1,0},{0,1},{1,1}};
    double g2_t[2][3] = {{1,0,1},{0,1,1}};
    double b2[3][3];
    double cov2[3][3];
    double cov2_inv[3][3];
    double mat_temp3[2][3];
    double mat_temp4[2][2];
    double mat_temp4_inv[2][2];
    double det;
    double ex_temp;
    double mat_temp5[2][3];
    double theta2[2] = {0,0};   

    //bs -> bs_mat;
    for (i=0; i<bs_row; i++){
        bs_mat[i][0] = *(bs+2*i);
        bs_mat[i][1] = *(bs+2*i+1);
    }

    //cal k1,g1,h1;
    for (i=0; i<bs_row; i++){
        k1[i] = pow(bs_mat[i][0],2) + pow(bs_mat[i][1],2);
        g1[i][0] = (-2)*bs_mat[i][0];
        g1[i][1] = (-2)*bs_mat[i][1];
        g1[i][2] = 1;
        h1[i] = pow(*(dis+i),2)-k1[i];
    }

    //g1 transpose;
    for (i=0; i<bs_row; i++){
        g1_t[0][i] = g1[i][0];
        g1_t[1][i] = g1[i][1];
        g1_t[2][i] = g1[i][2];
    }

    //g1_t*q1_inv;
    for (i=0; i<3; i++){
        for (j=0; j<bs_row; j++){
            mat_temp0[i][j] = g1_t[i][j]*(1/sigma);
        }
    }

    //mat_temp0*g1;
    for (i=0; i<3; i++){
        for (j=0; j<bs_row; j++){            
            mat_temp1[i][0] += mat_temp0[i][j] * g1[j][0];
            mat_temp1[i][1] += mat_temp0[i][j] * g1[j][1];
            mat_temp1[i][2] += mat_temp0[i][j] * g1[j][2];  
        }
    }
    for (i=0; i<3; i++){
        for (j=0; j<bs_row; j++){
            mat_temp0[i][j] = 0;
        }
    }

    //mat_temp1_inv;
    row_col = 3;
    arr=(double*)malloc(sizeof(double)*row_col*row_col);
    result=(double*)malloc(sizeof(double)*row_col*row_col);
    for (i=0; i<row_col; i++){
        for (j=0; j<row_col; j++){
            arr[i*row_col+j] = mat_temp1[i][j];
        }
    }
    MatrixOpp(arr, row_col, row_col, result);
    for (i=0; i<row_col; i++){
            for (j=0; j<row_col; j++){
                mat_temp1_inv[i][j] = *(result+i*3+j);
            }
    }
    free(result);
    free(arr);
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            mat_temp1[i][j] = 0;
        }
    }

    //mat_temp1_inv*g1_t;
    for (i=0; i<3; i++){
        for (j=0; j<bs_row; j++){
        	sum = 0;
            for (k=0; k<3; k++){
                sum += mat_temp1_inv[i][k] * g1_t[k][j];
            }
            mat_temp0[i][j] = sum;
        }
    }
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            mat_temp1_inv[i][j] = 0;
        }
    }

    //mat_temp0*q1_inv;
    for (i=0; i<3; i++){
        for (j=0; j<bs_row; j++){
            mat_temp0[i][j] = mat_temp0[i][j] * (1/sigma);
        }
    }

    //theta0: first WLS;
    for (i=0; i<bs_row; i++){
        theta0[0] += mat_temp0[0][i]*h1[i];
        theta0[1] += mat_temp0[1][i]*h1[i];
    }
    for (i=0; i<3; i++){
        for (j=0; j<bs_row; j++){
            mat_temp0[i][j] = 0;
        }
    }

    //cal the true distance;
    for (i=0; i<bs_row; i++){
        dis_t[i] = pow((bs_mat[i][0]-theta0[0]),2) + pow((bs_mat[i][1]-theta0[1]),2);
        dis_t[i] = sqrt(dis_t[i]);
    }

    //diag b1;
    for (i=0; i<bs_row; i++){
        for (j=0; j<bs_row; j++){
            b1[i][j] = 0;
        }
    }
    for (i=0; i<bs_row; i++){
        b1[i][i] = dis_t[i];
    }

    //b1*q1*b1 -> cov1;
    for (i=0; i<bs_row; i++){
        for (j=0; j<bs_row; j++){
            cov1[i][j] = 0;
        }
    }
    for (i=0; i<bs_row; i++){
        cov1[i][i] = b1[i][i]*sigma*b1[i][i];
    }

    //cov1_inv;
    for (i=0; i<bs_row; i++){
        for (j=0; j<bs_row; j++){
            cov1_inv[i][j] = 0;
        }
    }
    for (i=0; i<bs_row; i++){
        cov1_inv[i][i] = 1/cov1[i][i];
    }

    //g1_t*cov1_inv;
    for (i=0; i<3; i++){
        for (j=0; j<bs_row; j++){
        	sum = 0;
            for (k=0; k<bs_row; k++){
                sum += g1_t[i][k]*cov1_inv[k][j];
            }
            mat_temp0[i][j] = sum;
        }
    }

    //mat_temp0*g1;
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
        	sum = 0;
            for (k=0; k<bs_row; k++){
                sum += mat_temp0[i][k]*g1[k][j];
            }
            mat_temp1[i][j] = sum;
        }
    }
    for (i=0; i<3; i++){
        for (j=0; j<bs_row; j++){
            mat_temp0[i][j] = 0;
        }
    }

    //mat_temp1_inv;
    row_col = 3;
    arr=(double*)malloc(sizeof(double)*row_col*row_col);
    result=(double*)malloc(sizeof(double)*row_col*row_col);
    for (i=0; i<row_col; i++){
        for (j=0; j<row_col; j++){
            arr[i*row_col+j] = mat_temp1[i][j];
        }
    }
    MatrixOpp(arr, row_col, row_col, result);
    for (i=0; i<row_col; i++){
            for (j=0; j<row_col; j++){
                mat_temp1_inv[i][j] = *(result+i*3+j);
            }
    }
    free(result);
    free(arr);
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            mat_temp1[i][j] = 0;
        }
    }

    //mat_temp1_inv*g1_t;
    for (i=0; i<3; i++){
        for (j=0; j<bs_row; j++){
        	sum = 0;
            for (k=0; k<3; k++){
                sum += mat_temp1_inv[i][k]*g1_t[k][j];
            }
            mat_temp0[i][j] = sum;
        }
    }

    //mat_temp0*cov1_inv;
    for (i=0; i<3; i++){
        for (j=0; j<bs_row; j++){
        	sum = 0;
            for (k=0; k<bs_row; k++){
                sum += mat_temp0[i][k]*cov1_inv[k][j];
            }
            mat_temp2[i][j] = sum;
        }
    }
    for (i=0; i<3; i++){
        for (j=0; j<bs_row; j++){
            mat_temp0[i][j] = 0;
        }
    }

    //theta1: second WLS;
    for (i=0; i<3; i++){
        for (j=0; j<bs_row; j++){
            theta1[i] += mat_temp2[i][j]*h1[j];       
        }
    }
    for (i=0; i<3; i++){
        for (j=0; j<bs_row; j++){
            mat_temp2[i][j] = 0;
        }
    }

    //cal h2;
    h2[0] = pow(theta1[0],2);
    h2[1] = pow(theta1[1],2);
    h2[2] = theta1[2];
    
    //diag b2;
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            b2[i][j] = 0;
        }
    }
    b2[0][0] = theta1[0];
    b2[1][1] = theta1[1];
    b2[2][2] = 0.5;

    //b2*mat_temp1_inv;
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
        	sum = 0;
            for (k=0; k<3; k++){
                sum += b2[i][k]*mat_temp1_inv[k][j];
            }
            mat_temp1[i][j] = sum;
        }
    }

    //mat_temp1*b2;
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
        	sum = 0;
            for (k=0; k<3; k++){
                sum += mat_temp1[i][k]*b2[k][j];
            }
            cov2[i][j] = sum;
        }
    }
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            cov2[i][j] =cov2[i][j]*4;
        }
    }
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            mat_temp1[i][j] = 0;
        }
    }

    //cov2_inv;
    row_col = 3;
    arr=(double*)malloc(sizeof(double)*row_col*row_col);
    result=(double*)malloc(sizeof(double)*row_col*row_col);
    for (i=0; i<row_col; i++){
        for (j=0; j<row_col; j++){
            arr[i*row_col+j] = cov2[i][j];
        }
    }
    MatrixOpp(arr, row_col, row_col, result);
    for (i=0; i<row_col; i++){
            for (j=0; j<row_col; j++){
                cov2_inv[i][j] = *(result+i*3+j);
            }
    }
    free(result);
    free(arr);

    //g2_t*cov2_inv;
    for (i=0; i<2; i++){
        for (j=0; j<3; j++){
        	sum = 0;
            for (k=0; k<3; k++){
                sum += g2_t[i][k]*cov2_inv[k][j];
            }
            mat_temp3[i][j] = sum;
        }
    }

    //mat_temp3*g2;
    for (i=0; i<2; i++){
        for (j=0; j<2; j++){
        	sum = 0;
            for (k=0; k<3; k++){
                sum += mat_temp3[i][k]*g2[k][j];
            }
            mat_temp4[i][j] = sum;
        }
    }
    for (i=0; i<2; i++){
        for (j=0; j<3; j++){
            mat_temp3[i][j] = 0;
        }
    }

    //mat_temp4_inv;
    det = mat_temp4[0][0]*mat_temp4[1][1] - mat_temp4[0][1]*mat_temp4[1][0];
    ex_temp = mat_temp4[0][0];
    mat_temp4[0][0] = mat_temp4[1][1];
    mat_temp4[1][1] = ex_temp;
    mat_temp4[0][1] = -mat_temp4[0][1];
    mat_temp4[1][0] = -mat_temp4[1][0];
    for (i=0; i<2; i++){
        for (j=0; j<2; j++){
            mat_temp4_inv[i][j] = mat_temp4[i][j]*(1/det);
        }
    }
    for (i=0; i<2; i++){
        for (j=0; j<2; j++){
            mat_temp4[i][j] = 0;
        }
    }

    //mat_temp4_inv*g2_t;
    for (i=0; i<2; i++){
        for (j=0; j<3; j++){
        	sum = 0;
            for (k=0; k<2; k++){
                sum += mat_temp4_inv[i][k]*g2_t[k][j];
            }
            mat_temp3[i][j] = sum;
        }
    }

    //mat_temp3*cov2_inv;
    for (i=0; i<2; i++){
        for (j=0; j<3; j++){
        	sum = 0;
            for (k=0; k<3; k++){
                sum += mat_temp3[i][k]*cov2_inv[k][j];
            }
            mat_temp5[i][j] = sum;
        }
    }

    //theta2: third WLS;
    for (i = 0; i < 3; i++){
        theta2[0] = theta2[0] + mat_temp5[0][i]*h2[i];
        theta2[1] = theta2[1] + mat_temp5[1][i]*h2[i];
    }

    //get theta;
    theta2[0] = round(sqrt(fabs(theta2[0])));
    theta2[1] = round(sqrt(fabs(theta2[1])));
    if (theta1[0]<0){
        theta2[0] = -theta2[0];
    }
    if (theta1[1]<0){
        theta2[1] = -theta2[1];
    }
    theta[0] = theta2[0];
    theta[1] = theta2[1];
}

/*----------------------------------------------------------------------------------------*/
/*--------------------------------------MAT_INV-------------------------------------------*/
/*----------------------------------------------------------------------------------------*/
void MatrixOpp(double A[], int m, int n, double* invmat){
    int i, j, x, y, k;
    double *SP = NULL, *AB = NULL, *B = NULL, X;
    SP = (double*)malloc(m*n*sizeof(double));
    AB = (double*)malloc(m*n*sizeof(double));
    B = (double*)malloc(m*n*sizeof(double));
    X = Surplus(A, m, n);
    X = 1/X;
 
    for (i=0; i<m; i++){
        for (j=0; j<n; j++){
            for (k=0; k<m*n; k++){
                B[k] = A[k];
                for (x=0; x<n; x++)
                    B[i*n+x] = 0;
                for (y=0; y<m; y++)
                    B[m*y+j] = 0;
                B[i*n+j] = 1;
                SP[i*n+j] = Surplus(B, m, n);
                AB[i*n+j] = X*SP[i*n+j];
            }
        }
    }
    MatrixInver(AB, m, n, invmat);
    free(SP);
    free(AB);
    free(B);
}
 
void MatrixInver(double A[], int m, int n, double* invmat){
    int i, j;
    double *B = invmat; 
    for (i=0; i<n; i++)
        for (j=0; j<m; j++)
            B[i*m+j] = A[j*n+i];
}
 
double Surplus(double A[], int m, int n){
    int i, j, k, p, r;
    double X, temp=1, temp1=1, s=0, s1=0;
    if (n==2){
        for (i=0; i<m; i++)
            for (j=0; j<n; j++)
                if ((i+j)%2)
                    temp1 *= A[i*n+j];
                else
                    temp *= A[i*n+j];
        X = temp-temp1;
    }
    else{
        for (k=0; k<n; k++){
            for (i=0, j=k; i<m, j<n; i++, j++)
                temp *= A[i*n+j];
            if (m-i){
                for (p=m-i, r=m-1; p>0; p--, r--)
                    temp *= A[r*n+p-1];
            }
            s += temp;
            temp = 1;
        }
        for (k=n-1; k>=0; k--){
            for (i=0, j=k; i<m, j>=0; i++, j--)
                temp1 *= A[i*n+j];
            if (m-i){
                for (p=m-1, r=i; r<m; p--, r++)
                    temp1 *= A[r*n+p];
            }
            s1 += temp1;
            temp1 = 1;
        }
        X = s-s1;
    }
    return X;
}