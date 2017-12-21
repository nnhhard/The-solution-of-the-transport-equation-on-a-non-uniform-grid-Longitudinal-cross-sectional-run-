#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "Solvers.h"

const double a_x = 0, b_x = 1;
const double a_y = 0, b_y = 1;

const int N = 10;
const int M = 10;
const double T = 1000;

const double TAU = 0.01;
const double a = 1;

/*
double Derivative_x(double **vector,int i,int j,double *hx)
{
    return (vector[i+1][j] - vector[i-1][j]) / (hx[i-1] * hx[i]);
}
double Derivative_y(double **vector,int i,int j,double *hy)
{
    return (vector[i][j+1] - vector[i][j-1]) / (hy[i-1] * hy[i]);
}
//*/
double Laplas_x(double **vector,int i,int j,double *hx)
{
    return ( ( ( ( vector[i+1][j] - vector[i][j] ) / (hx[i]) ) - ( ( vector[i][j] - vector[i-1][j] ) / (hx[i-1] ) ) ) / ( ( hx[i] + hx[i-1] ) / 2.0 ) );
}
double Laplas_y(double **vector,int i,int j,double *hy)
{
    return ( ( ( ( vector[i][j+1] - vector[i][j] ) / (hy[j]) ) - ( ( vector[i][j] - vector[i][j-1] ) / (hy[j-1] ) ) ) / ( ( hy[j] + hy[j-1] ) / 2.0 ) );
}
void null(double **vector, int lN, int lM)
{
    for(int i = 0; i <= lN; i++)
        for(int j = 0; j <= lM; j++)
            vector[i][j]=0;
}
double Norm(double **vector,int lN,int lM)
{
    double max=fabs(vector[0][0]);

    for(int i=0;i<=lN;i++)
    {
        for(int j=0;j<=lM;j++)
        {
            if(fabs(vector[i][j])>max)
            {
                max=fabs(vector[i][j]);
            }
        }
    }
    return max;
}
double Norm(double *vector,int lN)
{
    double max = fabs(vector[0]);

    for(int i = 0; i <= lN; i++)
    {
        if(fabs(vector[i])>max)
        {
            max=fabs(vector[i]);
        }
    }
    return max;
}
double Function(double t, double x, double y)
{
    //return sin(M_PI * t) + sin(M_PI * x) * sin(M_PI * y);
    //return sin(M_PI * x) * sin(M_PI * y);
    //return t*t + sin(M_PI * x) * sin(M_PI * y);
    //return t * t + x * x + y * y;
    return x * x + y * y;
}
double RightPart(double t, double x, double y)
{
    //return cos(M_PI * t) + a * ( 2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y) );
    //return a * ( 2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y) );
    //return 2 * t + a * ( 2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y) );
    //return 2 * t - a * 4;
    return - a * 4;
}
void SolveTransport(double **u,double **u_n,double **mask,int lN,int lM,double *x,double *y,double *hx,double *hy,double t_1_2,double t)
{
    double *u_temp_x=new double[lN + 1];
    double *u_temp_y=new double[lM + 1];

    double **u_1_2=new double*[lN + 1];
    for(int i=0;i<=lN;i++)
        u_1_2[i]=new double[lM + 1];


    for(int i = 0; i <= lN; i++) {
        for(int j = 0; j <= lM; j++) {
            u_1_2[i][j] = u_n[i][j];
        }
    }

    double *A=new double[lM + 1];
    double *B=new double[lM + 1];
    double *C=new double[lM + 1];
    double *F=new double[lM + 1];

    for(int j = 0; j <= lM; j++)
    {
        A[j] = 0;
        B[j] = 0;
        C[j] = 0;
        F[j] = 0;
    }

    for(int i = 1; i <= lN - 1; i++)
    {
        for(int j = 0; j<= lM; j++)
        {
                if(j == 0)
                {
                    A[j]=0;
                    B[j]=1;
                    C[j]=0;
                    F[j]=Function(t_1_2,x[i],y[j]);
                }
                else if(j == lM)
                {
                    A[j]=0;
                    B[j]=1;
                    C[j]=0;
                    F[j]=Function(t_1_2,x[i],y[j]);
                }
                else
                {
                    A[j] =  ( - 1.0 / (hy[j-1]) ) / ( (hy[j] + hy[j-1]) / 2.0 ) * a;
                    B[j] =  2.0/TAU
                    + a
                    *(
                        (1.0 / hy[j] + 1.0 / hy[j-1]) / ( (hy[j]+hy[j-1]) / 2.0 )
                     );
                    C[j] =  ( - 1.0 / (hy[j]) ) / ( (hy[j]+hy[j-1]) / 2.0 ) * a;
                    F[j] = RightPart(t_1_2,x[i],y[j]) + a * Laplas_x(u_n,i,j,hx) + 2.0*u_n[i][j]/TAU;
                }
        }
        SolveByScalarRun(lM,u_temp_y,A,B,C,F);

        for(int j=0;j<=lM;j++)
            u_1_2[i][j]=u_temp_y[j];

        }

        for(int i=0;i<=lN;i++)
        {
            for(int j=0;j<=lM;j++)
            {
                u[i][j] = u_1_2[i][j];
            }
        }

        delete []A;
        delete []B;
        delete []C;
        delete []F;


        A=new double[lN + 1];
        B=new double[lN + 1];
        C=new double[lN + 1];
        F=new double[lN + 1];

        for(int j=1;j<lM;j++)
        {
            for(int i = 0; i<=lN; i++)
            {
                if(i==0)
                {
                    A[i]=0;
                    B[i]=1;
                    C[i]=0;
                    F[i]=Function(t,x[i],y[j]);
                }
                else if(i==lN)
                {
                    A[i]=0;
                    B[i]=1;
                    C[i]=0;
                    F[i]=Function(t,x[i],y[j]);
                }
                else
                {
                    A[i] =  ( - 1.0 / (hx[i-1]) ) / ( (hx[i] + hx[i-1]) / 2.0 ) * a;
                    B[i] =  2.0/TAU
                    + a
                    *(
                        (1.0 / hx[i] + 1.0 / hx[i-1]) / ( (hx[i] + hx[i-1]) / 2.0 )
                     );
                    C[i] =  ( - 1.0 / (hx[i]) ) / ( (hx[i] + hx[i-1]) / 2.0 ) * a;
                    F[i] = RightPart(t_1_2,x[i],y[j]) + a * Laplas_y(u_1_2,i,j,hy) + 2.0*u_1_2[i][j]/TAU;
                }
            }
        SolveByScalarRun(lN,u_temp_x,A,B,C,F);

        for(int i=0;i<=lN;i++)
            u[i][j]=u_temp_x[i];
        }

        for(int i=1;i<lN;i++)
        {
            u[i][0] = Function(t,x[i],y[0]);
            u[i][lM] = Function(t,x[i],y[lM]);
        }


    for(int i=0;i<=lN;i++)
    {
        delete []u_1_2[i];
    }
    delete []u_1_2;
    delete []A;
    delete []B;
    delete []C;
    delete []F;
    delete []u_temp_x;
    delete []u_temp_y;
    u_1_2=NULL;
    A=NULL;
    B=NULL;
    C=NULL;
    F=NULL;
    u_temp_x=NULL;
    u_temp_y=NULL;
}



int main()
{
    double *hx = new double[N];
    double *hy = new double[M];
    double h = (b_x - a_x) / (N);

    if(0 == 0)
    {
        hx[0] = h/4.0;
        hx[N-1] = h/4.0;

        for(int i = 1; i <= N-2; i++)
        {
            hx[i] = ( b_x - h/2.0 - a_x )/(N-2);
        }
    }
    else
    {
        for(int i = 0; i <= N-1; i++)
        {
            hx[i] = h;
        }
    }

    h = (b_y - a_y) / (M);
    for(int i = 0; i <= M-1; i++)
    {
        hy[i] = h;
    }

    double *x = new double[N+1];

    x[0]=a_x;
    for(int i = 1; i <= N; i++)
    {
        x[i] = x[i-1] + hx[i-1];
    }

    double *y = new double[M+1];

    y[0]=a_y;
    for(int i = 1; i <= M; i++)
    {
        y[i] = y[i-1] + hy[i-1];
    }


    double **u = new double*[N + 1];
    for(int i = 0; i <= N; i++)
        u[i] = new double[M + 1];

    double **mask_u = new double*[N + 1];
    for(int i = 0; i <= N; i++)
        mask_u[i] = new double[M + 1];

    double **u_n = new double*[N + 1];
    for(int i = 0; i <= N; i++)
        u_n[i] = new double[M + 1];

    double **w = new double*[N + 1];
    for(int i = 0; i <= N; i++)
        w[i] = new double[M + 1];


    null(w,N,M);
    null(u,N,M);
    null(mask_u,N,M);
    null(u_n,N,M);

     for(int i = 1; i <= N-1; i++){
        for(int j = 1; j <= M-1; j++){
            mask_u[i][j] = 1;
        }
    }

    for(int i = 0; i <= N; i++){
        for(int j = 0; j <= M; j++){
            u[i][j] = Function(0,x[i],y[j]);
            u_n[i][j] = Function(0,x[i],y[j]);
        }
    }

    for(int t=0;t<=T;t++)
    {
        //printf("%lf \n",((double)t + 1.0)*TAU);

        SolveTransport(u,u_n,mask_u,N,M,x,y,hx,hy,((double)t + 1.0/2.0)*TAU,((double)t + 1.0)*TAU);

        for(int i = 0; i <= N; i++) {
            for(int j = 0; j <= M; j++) {
                if ( (i==0) || (i==N) || (j==0) || (j==M) ) {
                    w[i][j]=0;
                }
                else {
                    w[i][j]=u[i][j]-Function(((double)t + 1.0)*TAU,x[i],y[j]);
                }
            }
        }
        printf("Norm(w) = %lf\n",Norm(w,N,M));
        for(int i = 0; i <= N; i++) {
            for(int j = 0; j <= M; j++) {
                u_n[i][j] = u[i][j];
            }
        }

    }

    for(int i = 0; i <= N; i++){
        delete []u[i];
        delete []u_n[i];
        delete []w[i];
        delete []mask_u[i];
    }
    delete []u;
    delete []u_n;
    delete []w;
    delete []mask_u;

    u = NULL;
    w = NULL;
    u_n = NULL;
    mask_u = NULL;

    return 0;
}

