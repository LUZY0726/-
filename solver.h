#ifndef __SOVLER_ONE_H__
#define __SOVLER_ONE_H__

#include "tri_block.h"
#include <math.h>

void solver_one(int scale)//纵向矩形块分块
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

void solver_two(int scale)//横向矩形块分块
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

void solver_three(int scale)//固定大小矩形块分块
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

void solver_four(int scale)//类似于递归的分块
{
    int ser_block=pow(2,scale);//串行三角的块数
    int size=m/ser_block;//每一三角块的行数
    if(size==0)
    size=1;//最小块大小为1

    int *figure;//设置2的次方数组方便rec_down的计算
    figure=(int *)malloc(sizeof(int)*scale+1);

    for(int i=0;i<=scale;i++)//初始化figure
    figure[i]=(int)pow(2,i);

    int *judge;//设置数组方便rec_down的计算
    judge=(int *)malloc(sizeof(int)*scale);
    memset(judge,0,sizeof(int)*scale);

    for(int i=0;i<scale;i++)//初始化
    {
        for(int j=scale-1;j>=scale-1-i;j--)
        {
            judge[i]+=figure[j];
        }
    }

    int tri_up=0;//初始化各参数值
    int tri_down=size;
    int rec_down;
    int rec_left;
    int rec_right;

    int count=0;//计算进行到了第几块

    while(count<ser_block-1)//只进行到ser_block-1块
    {
        // printf("tri_up=%d\n",tri_up);
        // printf("tri_down=%d\n",tri_down);

        part_serial(tri_up,tri_down);//先进行三角部分串行运算
        count++;//目前进行的块数

        int flag=0;//标志数
        int temp=count;//方便将count的值进行计算
        while(1)
        {
            for(int i=0;i<=scale;i++)
            {
                if(temp==figure[i])//如果count值等于2的次方则计算出矩形的下界值
                {
                    rec_down=tri_down+temp*size;
                    flag=1;
                    break;
                }
            }

            if(flag==1)//下界值求出死循环结束
            break;

            for(int i=0;i<=scale;i++)//如果count值不等于2的次方则做减法运算
            {
                if(temp<figure[i])
                {
                    temp-=figure[i-1];
                    break;
                }
            }
        }

        rec_right=tri_down;//矩形块的右界值等于对应三角块的下界值
        rec_left=rec_right-(rec_down-tri_down);//矩形的宽度等于深度
        // printf("rec_left=%d\n",rec_left);
        // printf("rec_right=%d\n",rec_right);

        for(int i=0;i<scale;i++)
        {
            if(judge[i]==count)//如果count值等于judge数组中的任意一个则下界为矩阵底
            {
                rec_down=m;
                break;
            }
        }
        //printf("rec_down=%d\n\n",rec_down);

        part_parallel(tri_down,rec_down,rec_left,rec_right);//求解矩形块

        tri_up+=size;//三角部分参数自加
        tri_down+=size;
    }
    part_serial(tri_up,m);//最后一个三角块单独运算且无此时已无矩形块
    // printf("tri_up=%d\n",tri_up);
    // printf("tri_down=%d\n",m);
}

#endif
