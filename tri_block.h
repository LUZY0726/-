#ifndef __TRI_BLOCK_H__
#define __TRI_BLOCK_H__

#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include "mmio_highlevel.h"

int *csrRowPtrA;  //原矩阵CSR行信息
int *csrColIdxA;  //原矩阵CSR列下标
double *csrValA;   //原矩阵CSR非零元值
int *csrRowPtr_tmp;	 //三角解CSR行信息
int *csrColIdx_tmp;  //三角解CSR列下标
double *csrVal_tmp;  //三角解CSR非零元值

double *Result;      //结果向量
double *X;          //待求解向量
int m;     //矩阵行数
int n;		//矩阵列数
int nnzA;	//矩阵非零元数
int isSymmetricA;     //m是否等于n(矩阵是否对称)

// "Usage: ``./sptrsv -d 0 -rhs 1 -forward -mtx A.mtx'' for LX=B on device 0"
int argi = 1;
int substitution = SUBSTITUTION_FORWARD;
int nnzTR;

void part_serial(int ubound,int lbound)
{
	for(int i=ubound;i<lbound;i++)//将三角矩阵的上半部分作为第一部分优先求解 
	{
		double sumf=0;
		for(int j=csrRowPtr_tmp[i+1]-2;j>=csrRowPtr_tmp[i];j--)//遍历非零元素 
		{
			if(csrColIdx_tmp[j]>=ubound)
			sumf+=csrVal_tmp[j]*X[csrColIdx_tmp[j]];
			else
			{
				break;
			}
		}
			X[i]=(Result[i]-sumf)/csrVal_tmp[csrRowPtr_tmp[i+1]-1];//求解得出对应的目标向量元素
	}
 }

 void part_parallel(int ubound,int lbound,int left,int right)
 {
	#pragma omp parallel for
 	for(int i=ubound;i<lbound;i++)
	 {
		 for(int j=csrRowPtr_tmp[i];j<csrRowPtr_tmp[i+1];j++)
		 {
				 if(csrColIdx_tmp[j]<right&&csrColIdx_tmp[j]>=left)
				{
					Result[i]-= csrVal_tmp[j]*X[csrColIdx_tmp[j]];//结果向量减去矩形内的非零元
			 	}
		 }
	 }	 
 }

 #endif 

 
