#include "solver.h"
 
int main(int argc, char ** argv)
{
    char  *filename;
    if(argc > argi)
    {
        filename = argv[argi];
        argi++;
    }
    printf("\n-------------- %s --------------\n", filename);
	// load mtx data to the csr format
        mmio_info(&m, &n, &nnzA, &isSymmetricA, filename);
	//	printf("m=%d,n=%d,nnzA=%d\n", m, n, nnzA);
        csrRowPtrA = (int *)malloc((m+1) * sizeof(int));
        csrColIdxA = (int *)malloc(nnzA * sizeof(int));
        csrValA    = (VALUE_TYPE *)malloc(nnzA * sizeof(VALUE_TYPE));
        mmio_data(csrRowPtrA, csrColIdxA, csrValA, filename);
        printf("input matrix A: ( %i, %i ) nnz = %i\n", m, n, nnzA);
	//*******************************文件读取结束****************************************
	//*******************************将矩阵转化为L/U型***********************************
	// extract L or U with a unit diagonal of A
        csrRowPtr_tmp = (int *)malloc((m+1) * sizeof(int));
        csrColIdx_tmp = (int *)malloc((m+nnzA) * sizeof(int));
        csrVal_tmp    = (VALUE_TYPE *)malloc((m+nnzA) * sizeof(VALUE_TYPE));

        int nnz_pointer = 0;
        csrRowPtr_tmp[0] = 0;
        for (int i = 0; i < m; i++)
        {
            for (int j = csrRowPtrA[i]; j < csrRowPtrA[i+1]; j++)
            {   
                if (substitution == SUBSTITUTION_FORWARD)
                {
                    if (csrColIdxA[j] < i)
                    {
                        csrColIdx_tmp[nnz_pointer] = csrColIdxA[j];
                        csrVal_tmp[nnz_pointer] = 1;//rand() % 10 + 1; //csrValA[j]; 
                        nnz_pointer++;
                    }
                }
                else if (substitution == SUBSTITUTION_BACKWARD)
                {
                    if (csrColIdxA[j] > i)
                    {
                        csrColIdx_tmp[nnz_pointer] = csrColIdxA[j];
                        csrVal_tmp[nnz_pointer] = 1;//rand() % 10 + 1; //csrValA[j]; 
                        nnz_pointer++;
                    }
                }
            }

            // add dia nonzero
            csrColIdx_tmp[nnz_pointer] = i;
            csrVal_tmp[nnz_pointer] = 1.0;
            nnz_pointer++;

            csrRowPtr_tmp[i+1] = nnz_pointer;
        }
		//**************************************************8
		 int nnz_tmp = csrRowPtr_tmp[m];
        nnzTR = nnz_tmp;

        if (substitution == SUBSTITUTION_FORWARD)
            printf("A's unit-lower triangular L: ( %i, %i ) nnz = %i\n", m, n, nnzTR);
        else if (substitution == SUBSTITUTION_BACKWARD)
            printf("A's unit-upper triangular U: ( %i, %i ) nnz = %i\n", m, n, nnzTR);

    //  csrColIdx_tmp = (int *)realloc(csrColIdx_tmp, sizeof(int) * nnzTR);
    //  csrVal_tmp = (VALUE_TYPE *)realloc(csrVal_tmp, sizeof(VALUE_TYPE) * nnzTR);
	
	//*************************************转化完成************************************
	free(csrRowPtrA);
    free(csrColIdxA);
    free(csrValA);

    int scale=10;//设置三角块数pow(2,scale)个
    int flag=1;//计算正确标志位

    Result=(double *)malloc(sizeof(double)*m);
	memset(Result,0,sizeof(double)*m);

	//*****************************计算结果向量方便验证*****************************
	for(int i=0;i<m;i++)
	{
		for(int j=csrRowPtr_tmp[i];j<csrRowPtr_tmp[i+1];j++)
		{
			Result[i]+=csrVal_tmp[j];
		}
	}
	//***************************************************************************
	X=(double *)malloc(sizeof(double)*n);
    solver_three(scale);
	for(int i=0;i<n;i++)
	{
		if(X[i]!=1)
		{
			flag=0;
		}
	}
    if(flag==0)
    printf("error!!!\n");
    else
    printf("success!!!\n");

    // for(int i=0;i<n;i++)
    // printf("%.1lf ",X[i]);

    free(Result);
    free(X);
    free(csrRowPtr_tmp);
	free(csrColIdx_tmp);
	free(csrVal_tmp);
    return 0;
}
