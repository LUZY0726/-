#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <malloc.h>
#include "mmio_highlevel.h"

void printmat(float *A, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
            printf("%4.2f ", A[i * n + j]);
        printf("\n");
    }
}

void printvec(float *x, int n)
{
    for (int i = 0; i < n; i++)
        printf("%4.2f\n", x[i]);
}

void matvec(float *A, float *x, float *y, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        y[i] = 0;
        for (int j = 0; j < n; j++)
            y[i] += A[i * n + j] * x[j];
    }
}

void matmat(float *C, float *A, float *B, int m, int k, int n)
{
    memset(C, 0, sizeof(float) * m * n);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            for (int kk = 0; kk < k; kk++)
                C[i * n + j] += A[i * k + kk] * B[kk * n + j];
}

void matmat_transB(float *C, float *A, float *BT, int m, int k, int n)
{
    memset(C, 0, sizeof(float) * m * n);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            for (int kk = 0; kk < k; kk++)
                C[i * n + j] += A[i * k + kk] * BT[j * k + kk];
}
int count(float *A,float *BT,int m,int k,int n)
{
    int nnz=0;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            float sum=0;
            for (int kk = 0; kk < k; kk++)
            {
                sum+=A[i * k + kk] * BT[j * k + kk];
            }
            if(sum!=0)
            {
                nnz++;
            }
        }        
    }
    return nnz;
}
void matmat_transB_temp(int *csrrowptrrt,int *csrcolidxrt,float *csrvalrt,float *A,float *BT,int m,int k,int n)
{
    int count=0;
    csrrowptrrt[0]=0;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            float sum=0;
            for (int kk = 0; kk < k; kk++)
            {
                sum+=A[i * k + kk] * BT[j * k + kk];
            }
            if(sum!=0)
            {
                csrcolidxrt[count]=j;
                csrvalrt[count]=sum;
                count++;
            }
        }
        csrrowptrrt[i+1]=count;
    }
}
float dotproduct(float *vec1, float *vec2, int n)
{
    float result = 0;
    for (int i = 0; i < n; i++)
        result += vec1[i] * vec2[i];
    return result;
}

// A is m x n, AT is n x m
void transpose(float *AT, float *A, int m, int n)
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            AT[j * m + i] = A[i * n + j];
}

float vec2norm(float *x, int n)
{
    float sum = 0;
    for (int i = 0; i < n; i++)
        sum += x[i] * x[i];
    return sqrt(sum);
}

void cg(float *A, float *x, float *b, int n, int *iter, int maxiter, float threshold)
{
    memset(x, 0, sizeof(float) * n);
    float *residual = (float *)malloc(sizeof(float) * n);
    float *y = (float *)malloc(sizeof(float) * n);
    float *p = (float *)malloc(sizeof(float) * n);
    float *q = (float *)malloc(sizeof(float) * n);
    *iter = 0;
    float norm = 0;
    float rho = 0;
    float rho_1 = 0;

    // p0 = r0 = b - Ax0
    matvec(A, x, y, n, n);
    for (int i = 0; i < n; i++)
        residual[i] = b[i] - y[i];
    //printvec(residual, n);

    do
    {
        //printf("\ncg iter = %i\n", *iter);
        rho = dotproduct(residual, residual, n);
        if (*iter == 0)
        {
            for (int i = 0; i < n; i++)
                p[i] = residual[i];
        }
        else
        {
            float beta = rho / rho_1;
            for (int i = 0; i < n; i++)
                p[i] = residual[i] + beta * p[i];
        }

        matvec(A, p, q, n, n);
        float alpha = rho / dotproduct(p, q, n);
        //printf("alpha = %f\n", alpha);
        for (int i = 0; i < n; i++)
            x[i] += alpha * p[i];
        for (int i = 0; i < n; i++)
            residual[i] += - alpha * q[i];

        rho_1 = rho;
        float error = vec2norm(residual, n) / vec2norm(b, n);

        //printvec(x, n);
        *iter += 1;

        if (error < threshold)
            break;
    }
    while (*iter < maxiter);

    free(residual);
    free(y);
    free(p);
    free(q);
}

void updateX(int *csrRowPtrR, int *csrColIdxR, float *csrValR, float *X, float *Y,
             int m, int n, int f, float lamda)
{
    // create YT, and transpose Y
    float *YT = (float *)malloc(sizeof(float) * n * f);
    transpose(YT, Y, n, f);

    // create smat = YT*Y
    float *smat = (float *)malloc(sizeof(float) * f * f);

    // multiply YT and Y to smat
    matmat(smat, YT, Y, f, n, f);

    // smat plus lamda*I
    for (int i = 0; i < f; i++)
        smat[i * f + i] += lamda;

    // create svec
    float *svec = (float *)malloc(sizeof(float) * f);

    // loop for all rows of X, from 0 to m-1
    for (int u = 0; u < m; u++)
    {
        printf("\n u = %i", u);
        // reference the uth row of X
        float *xu = &X[u * f];

        // compute svec by multiplying YT and the uth row of R
        float *temp=(float *)malloc(sizeof(float)*n);
        memset(temp,0,sizeof(float)*n);
        for(int i=csrRowPtrR[u];i<csrRowPtrR[u+1];i++)
        {
            temp[csrColIdxR[i]]=csrValR[i];
        }
        matvec(YT, temp, svec, f, n);

        // solve the Ax=b system (A is smat, b is svec, x is xu (uth row of X))
        int cgiter = 0;
        cg(smat, xu, svec, f, &cgiter, 100, 0.00001);

        free(temp);
        //printf("\nsmat = \n");
        //printmat(smat, f, f);

        //printf("\nsvec = \n");
        //printvec(svec, f);

        //printf("\nxu = \n");
        //printvec(xu, f);
    }

    free(smat);
    free(svec);
    free(YT);
}

void updateY(int *csrRowPtrR, int *csrColIdxR, float *csrValR, float *X, float *Y,
             int m, int n, int f, float lamda)
{
    float *XT = (float *)malloc(sizeof(float) * m * f);
    transpose(XT, X, m, f);

    float *smat = (float *)malloc(sizeof(float) * f * f);
    matmat(smat, XT, X, f, m, f);
    for (int i = 0; i < f; i++)
        smat[i * f + i] += lamda;

    float *svec = (float *)malloc(sizeof(float) * f);
    float *ri = (float *)malloc(sizeof(float) * m);

    for (int i = 0; i < n; i++)
    {
        printf("\n i = %i", i);
        float *yi = &Y[i * f];

        memset(ri,0,sizeof(float)*m);
        for (int k = 0; k < m; k++)
        {
            //ri[k] = R[k * n + i];
            for(int j=csrRowPtrR[k];j<csrRowPtrR[k+1];j++)
            {
                if(csrColIdxR[j]==i)
                ri[k]=csrValR[j];
            }
        }
            

        matvec(XT, ri, svec, f, m);

        int cgiter = 0;
        cg(smat, yi, svec, f, &cgiter, 100, 0.00001);

        printf("\nsmat = \n");
        printmat(smat, f, f);

        printf("\nsvec = \n");
        printvec(svec, f);

        printf("\nyi = \n");
        printvec(yi, f);
    }

    free(smat);
    free(svec);
    free(ri);
    free(XT);
}

// ALS for matrix factorization
/*void als(int *csrRowPtrR, int *csrColIdxR, float *csrValR, float *X, float *Y,
         int m, int n, int f, float lamda)
{
    // create YT
    float *YT = (float *)malloc(sizeof(float) * f * n);

    // create R (result)
    float *Rp = (float *)malloc(sizeof(float) * m * n);

    int iter = 0;
    float error = 0.0;
    float error_old = 0.0;
    float error_new = 0.0;
    do
    {
        // step 1. update X
        updateX(csrRowPtrR, csrColIdxR, csrValR, X, Y, m, n, f, lamda);

        // step 2. update Y
        updateY(csrRowPtrR, csrColIdxR, csrValR, X, Y, m, n, f, lamda);

        // step 3. validate R, by multiplying X and YT

        // step 3-1. matrix multiplication with transposed Y
        matmat_transB(Rp, X, Y, m, f, n);

        // step 3-2. calculate error
        error_new = 0.0;
        for (int i = 0; i < m * n; i++)
            if (R[i] != 0) // only compare nonzero entries in R
                error_new += fabs(Rp[i] - R[i]) * fabs(Rp[i] - R[i]);
        error_new = sqrt(error_new/(m * n));

        error = fabs(error_new - error_old) / error_new;
        error_old = error_new;
        printf("iter = %i, error = %f\n", iter, error);

        iter++;
    }
    while(iter < 10 && error > 0.0001);

    printf("\nR = \n");
    printmat(R, m, n);

    printf("\nRp = \n");
    printmat(Rp, m, n);

    free(Rp);
    free(YT);
}*/

void updateX_recsys(int *csrRowPtrR, int *csrColIdxR, float *csrValR, float *X, float *Y,
                    int m, int n, int f, float lamda,
                    double *time_prepareA, double *time_prepareb, double *time_solver)
{
    struct timeval t1, t2;

    // malloc smat (A) and svec (b)
    float *smat = (float *)malloc(sizeof(float) * f * f);
    float *svec = (float *)malloc(sizeof(float) * f);

    for (int u = 0; u < m; u++)
    {
        gettimeofday(&t1, NULL);
        //printf("\n u = %i", u);
        float *xu = &X[u * f];

        // find nzr (i.e., #nonzeros in the uth row of R)
        int nzr = 0;
        // for (int k = 0; k < n; k++)
        //     nzr = R[u * n + k] == 0 ? nzr : nzr + 1;
        
        nzr = csrRowPtrR[u+1]-csrRowPtrR[u];//revise
        //printf("nzr=%d\n",nzr);

        // malloc ru (i.e., uth row of R) and insert entries into it
        float *ru = (float *)malloc(sizeof(float) * nzr);
        // int count = 0;
        // for (int k = 0; k < n; k++)
        // {
        //     if (R[u * n + k] != 0)
        //     {
        //         ru[count] = R[u * n + k];
        //         count++;
        //     }
        // }

        //revise
        for(int k = 0; k < nzr; k++)
        {
            ru[k] = csrValR[csrRowPtrR[u]+k];
        }
    /*    if(u<10)
        {
            for(int i=0;i<nzr;i++)
            printf("ru[%d]=%f ",i,ru[i]);
            printf("\n");
        }*/
        
        //printf("\n nzr = %i, ru = \n", nzr);
        //printvec(ru, nzr);

        // create sY and sYT (i.e., the zero-free version of Y and YT)
        float *sY = (float *)malloc(sizeof(float) * nzr * f);
        float *sYT = (float *)malloc(sizeof(float) * nzr * f);

        // fill sY, according to the sparsity of the uth row of R
        //int count = 0;
        // for (int k = 0; k < n; k++)
        // {
        //     if (R[u * n + k] != 0)
        //     {
        //         memcpy(&sY[count * f], &Y[k * f], sizeof(float) * f);
        //         count++;
        //     }
        // }

        //revise
        int count = 0;
        for(int k = 0; k < nzr; k++)
        {
            memcpy(&sY[count * f], &Y[csrColIdxR[csrRowPtrR[u]+k] * f], sizeof(float) * f);
            count++;
        //    if(u<10)
        //    printf("%d\n",csrColIdxR[csrRowPtrR[u]+k] *f);
        }
        //printf("\n sY = \n");
        //printmat(sY, nzr, f);

        // transpose sY to sYT
        transpose(sYT, sY, nzr, f);

        // multiply sYT and sY, and plus lamda * I
        matmat(smat, sYT, sY, f, nzr, f);
        for (int i = 0; i < f; i++)
            smat[i * f + i] += lamda;

        gettimeofday(&t2, NULL);
        *time_prepareA += (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;

        // compute b (i.e., svec) by multiplying sYT and the uth row of R
        gettimeofday(&t1, NULL);
        matvec(sYT, ru, svec, f, nzr);
        gettimeofday(&t2, NULL);
        *time_prepareb += (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;

        // solve the system of Ax=b, and get x = the uth row of X
        gettimeofday(&t1, NULL);
        int cgiter = 0;
        cg(smat, xu, svec, f, &cgiter, 100, 0.00001);
        gettimeofday(&t2, NULL);
        *time_solver += (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;

        //printf("\nsmat = \n");
        //printmat(smat, f, f);

        //printf("\nsvec = \n");
        //printvec(svec, f);

        //printf("\nxu = \n");
        //printvec(xu, f);

        free(ru);
        free(sY);
        free(sYT);
    }

    free(smat);
    free(svec);
    //free(YT);
}

void updateY_recsys(int *rowptrR,int *colidxR,float *valR,float *X,float *Y,int m,int n,int f,
float lamda,double *time_prepareA,double *time_prepareb,double *time_solver)
{
    struct timeval t1, t2;
    float *smat = (float *)malloc(sizeof(float) * f * f);
    float *svec = (float *)malloc(sizeof(float) * f);
        for(int i=0;i<n;i++)
        {
            gettimeofday(&t1, NULL);
            float *yi = &Y[i * f];
            int nzc=0;
            for(int k=0;k<rowptrR[m];k++)
            {
                if(colidxR[k]==i)
                nzc=nzc+1;
            }
        //    printf("nzc=%d\n",nzc);
            float *ri=(float *)malloc(sizeof(float)*nzc);
            int count=0;
            for(int l=0;l<rowptrR[m];l++)
            {
                if(colidxR[l]==i)
                {
                    ri[count]=valR[l];
                    count++;
                }
            }
        /*        if(i<10)
            {
                for(int k=0;k<count;k++)
                printf("ri[%d]=%f\n",k,ri[k]);
                printf("\n\n");
            }*/
            
    
    float *sX = (float *)malloc(sizeof(float) * nzc * f);
    float *sXT = (float *)malloc(sizeof(float) * nzc * f);
    count = 0;
         for(int l=0;l<rowptrR[m];l++)
         {
             if(colidxR[l]==i)
             {
                int k=0;
                for(int j=0;j<=m;j++)
                {
                    if(l>=rowptrR[j])
                    {
                        k=j;
                    }
                    else
                    {
                        break;
                    }
                    
                }
                memcpy(&sX[count * f], &X[k * f], sizeof(float) * f);
                count++;
            //    if(i<10)
            //    printf("%d\n",k * f);
             }
         }
    //     printf("\n\n");

    transpose(sXT, sX, nzc, f);
    matmat(smat, sXT, sX, f, nzc, f);
    for (int i = 0; i < f; i++)
            smat[i * f + i] += lamda;

        gettimeofday(&t2, NULL);
        *time_prepareA += (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;

        gettimeofday(&t1, NULL);
        matvec(sXT, ri, svec, f, nzc);
        gettimeofday(&t2, NULL);
        *time_prepareb += (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;

        gettimeofday(&t1, NULL);
        int cgiter = 0;

        

        cg(smat, yi, svec, f, &cgiter, 100, 0.00001);
     /*   if(i==1)
        {
            for(int h=0;h<f*f;h++)
            printf("%.1f ",smat[h]);
        }*/
        gettimeofday(&t2, NULL);
        *time_solver += (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;

        //printf("\nsmat = \n");
        //printmat(smat, f, f);

        //printf("\nsvec = \n");
        //printvec(svec, f);

        //printf("\nyi = \n");
        //printvec(yi, f);

        free(ri);
        free(sX);
        free(sXT);
        }
    free(smat);
    free(svec);
}

void als_recsys(int *rowptrR,int *colidxR,float *valR,float *X,float *Y,int m,int n,int f,float lamda)
{
    float *YT = (float *)malloc(sizeof(float) * f * n);
//    float *Rp = (float *)malloc(sizeof(float) * m * n);

    int iter = 0;
    float error = 0.0;
    float error_old = 0.0;
    float error_new = 0.0;
    struct timeval t1, t2;

    double time_updatex_prepareA = 0;
    double time_updatex_prepareb = 0;
    double time_updatex_solver = 0;

    double time_updatey_prepareA = 0;
    double time_updatey_prepareb = 0;
    double time_updatey_solver = 0;

    double time_updatex = 0;
    double time_updatey = 0;
    double time_validate = 0;
    do
    {
         // step 1. update X
        gettimeofday(&t1, NULL);
        updateX_recsys(rowptrR,colidxR,valR, X, Y, m, n, f, lamda,
                       &time_updatex_prepareA, &time_updatex_prepareb, &time_updatex_solver);
        gettimeofday(&t2, NULL);
        time_updatex += (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;

        // step 2. update Y
        gettimeofday(&t1, NULL);
        updateY_recsys(rowptrR,colidxR,valR, X, Y, m, n, f, lamda,
                       &time_updatey_prepareA, &time_updatey_prepareb, &time_updatey_solver);
        gettimeofday(&t2, NULL);
        time_updatey += (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;

        // step 3. validate
        // step 3-1. matrix multiplication
        int nnz1=count(X,Y,m,f,n);
        int *csrRowPtrRT=(int *)malloc(sizeof(int)*(m+1));
        int *csrColIdxRT=(int *)malloc(sizeof(int )*nnz1);
        float *csrValRT=(float *)malloc(sizeof(float)*nnz1);
        gettimeofday(&t1, NULL);
     //   matmat_transB(Rp, X, Y, m, f, n);
        
        matmat_transB_temp(csrRowPtrRT, csrColIdxRT, csrValRT,X, Y, m, f, n);
    //    printf("nnz1=%d\n",nnz1);
        // for(int i=0;i<nnz1;i++)
        // printf("%.1f ",csrValRT[i]);

        // step 3-2. calculate error
        error_new = 0.0;
        int nnz = 0;
//*********************************************************
        for(int i=0;i<m;i++)
        {
            for(int j=rowptrR[i];j<rowptrR[i+1];j++)
            {
                int count=0;
                for(int p=csrRowPtrRT[i];p<csrRowPtrRT[i+1];p++)
                {
                    count++;
                    if(csrColIdxRT[p]==colidxR[j])
                    {
                        error_new+=fabs(csrValRT[p]-valR[j])*fabs(csrValRT[p]-valR[j]);
                        break;
                    }
                    
                    if(count==csrRowPtrRT[i+1]-csrRowPtrRT[i])
                    {
                        error_new+=fabs(valR[j])*fabs(valR[j]);
                        break;
                    }   
                }      
                nnz++;
            }
        }
    /*    for(int i=0;i<rowptrR[m];i++)
        {
        //    printf("%f %f\n",csrValRT[i],valR[i]);
            
            error_new+=fabs(csrValRT[i]-valR[i])*fabs(csrValRT[i]-valR[i]);
            nnz++;
        }*/
        
        error_new = sqrt(error_new/nnz);
        gettimeofday(&t2, NULL);
        time_validate += (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;

        error = fabs(error_new - error_old) / error_new;
        error_old = error_new;
        printf("iter = %i, error = %f\n", iter, error);

        iter++;

        free(csrRowPtrRT);
        free(csrColIdxRT);
        free(csrValRT);

    } while (iter<1000 && error>0.0001);
     printf("\nUpdate X %4.2f ms (prepare A %4.2f ms, prepare b %4.2f ms, solver %4.2f ms)\n",
           time_updatex, time_updatex_prepareA, time_updatex_prepareb, time_updatex_solver);
    printf("Update Y %4.2f ms (prepare A %4.2f ms, prepare b %4.2f ms, solver %4.2f ms)\n",
           time_updatey, time_updatey_prepareA, time_updatey_prepareb, time_updatey_solver);
    printf("Validate %4.2f ms\n", time_validate);

   // free(Rp);
    
    free(YT);
    
}

int main(int argc, char ** argv)
{
    // parameters
    int f = 0;

    // method
    char *method = argv[1];
    printf("\n");

    char *filename = argv[2];
    printf ("filename = %s\n", filename);

    int m, n, nnzR, isSymmetricR;

    mmio_info(&m, &n, &nnzR, &isSymmetricR, filename);
    int *csrRowPtrR = (int *)malloc((m+1) * sizeof(int));
    int *csrColIdxR = (int *)malloc(nnzR * sizeof(int));
    float *csrValR    = (float *)malloc(nnzR * sizeof(float));
    mmio_data(csrRowPtrR, csrColIdxR, csrValR, filename);
//    for(int i=0;i<nnzR;i++)
//    printf("%f ",csrValR[i]);
    //FILE *file;
    //file = fopen(filename, "r+");
    //fscanf(file, "%i", &m);
    //fscanf(file, "%i", &n);
    printf("The order of the rating matrix R is %i by %i, #nonzeros = %i\n",
           m, n, nnzR);

    // create R
/*    float *R = (float *)malloc(sizeof(float) * m * n);
    memset(R, 0, sizeof(float) * m * n);*/

    // put nonzeros into the dense form of R (the array R)
/*    for (int rowidx = 0; rowidx < m; rowidx++)
    {
        for (int j = csrRowPtrR[rowidx]; j < csrRowPtrR[rowidx + 1]; j++)
        {
            int colidx = csrColIdxR[j];
            float val = csrValR[j];
            R[rowidx * n + colidx] = val;
        }
    }*/

    // read R
    //for (int i = 0; i < m; i++)
    //    for (int j = 0; j < n; j++)
    //        fscanf(file, "%f", &R[i * n + j]);

    //printf("\nR = \n");
    //printmat(R, m, n);

    f = atoi(argv[3]);
    printf("The latent feature is %i \n", f);

    // create X
    float *X = (float *)malloc(sizeof(float) * m * f);
    memset(X, 0, sizeof(float) * m * f);

    // create Y
    float *Y = (float *)malloc(sizeof(float) * n * f);
    for (int i = 0; i < n * f; i++)
        Y[i] = i%10;

    // lamda parameter
    float lamda = 0.1;

    // call function
/*    if (strcmp(method, "als") == 0)
    {
        als(csrRowPtrR,csrColIdxR,csrValR, X, Y, m, n, f, lamda);
    }*/
    if (strcmp(method, "als_recsys") == 0)
    {
        als_recsys(csrRowPtrR,csrColIdxR,csrValR, X, Y, m, n, f, lamda);
    }

 //   free(R);
    free(X);
    free(Y);

    return 0;
}

