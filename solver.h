#ifndef __SOVLER_ONE_H__
#define __SOVLER_ONE_H__

#include "tri_block.h"
#include <math.h>

void solver_one(int scale)//纵向
{
    int ser_block=pow(2,scale);//串行三角的块数
    int size=m/ser_block;//每一三角块的行数
    if(size==0)
    size=1;//最小块大小为1
    int count=0;//计算进行到多少块

    for(int i=0;i<m;i+=size)
    {
        count++;
        if(count==ser_block)//进行到最后一块的时候只进行串行解法并将多余行包含
        {
            part_serial(i,m);
            break;
        }
        else
        {
            part_serial(i,i+size);
            part_parallel(i+size,m,i,i+size);//左值left为上界i
        }    
    }
}

void solver_two(int scale)//横向
{
    int ser_block=pow(2,scale);//串行三角的块数
    int size=m/ser_block;//每一三角块的行数
    if(size==0)
    size=1;//最小块大小为1

    part_serial(0,size);//先进行第一块三角部分的求解
    int count=1;//计算进行到多少块
    for(int i=size;i<m;i+=size)
    {
        count++;
        if(count==ser_block)////进行到最后一部分的时候将多余行包含
        {
            part_parallel(i,m,0,i);
            part_serial(i,m);
            break;
        }
        part_parallel(i,i+size,0,i);
        part_serial(i,i+size);    
    }
}

void solver_three(int scale)//横向
{
    int ser_block=pow(2,scale);//串行三角的块数
    int size=m/ser_block;//每一三角块的行数
    if(size==0)
    size=1;//最小块大小为1

    part_serial(0,size);//先进行第一块三角部分的求解
    int count=1;//计算进行到多少块
    int block_num=1;//每一行的矩形块数目
    int left=0;//矩形块的左值
    for(int i=size;i<m;i+=size)
    {
        count++;
        if(count==ser_block)//进行到最后一部分的时候将多余行包含
        {
            left=0;
            for(int j=0;j<block_num;j++)
            {
                part_parallel(i,m,left,left+size);//每一次计算size列的矩形块
                left+=size;//左值变化
            }
            part_serial(i,m);
            break;
        }

        left=0;
        for(int j=0;j<block_num;j++)
        {
            part_parallel(i,i+size,left,left+size);//每一次计算size列的矩形块
            left+=size;//左值变化
        }
        part_serial(i,i+size);  

        block_num++;//每到下一行矩形块数目加1
    }
}

#endif
