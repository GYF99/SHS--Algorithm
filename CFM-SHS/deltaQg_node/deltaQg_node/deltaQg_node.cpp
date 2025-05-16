#include "D:\matlab\extern\include\mex.h"
#include <stdlib.h>
#include<string.h>
#include<stdio.h>
#include <ctype.h>
#include <vector>
#include <algorithm>   
#include <minmax.h>
#include <stdlib.h>
#include <vector>


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	//nlhs：输出参数个数
	//plhs：输出参数列表
	//nrhs：输入参数个数
	//prhs：输入参数列表   

	// Check for proper number of input and output arguments 
	if (nrhs != 8) {
		mexErrMsgIdAndTxt("MATLAB:mxcreatecellmatrix:nrhs",
			"输入参数个数不正确");
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("MATLAB:mxcreatecellmatrix:nlhs",
			"输出参数个数不正确");
	}

	//第一个输入变量的传递（原个体Ubefore）
	size_t mrowsU1;    //行数
	size_t ncolsU1;    //列数
	double *inMatU1;
	mrowsU1 = mxGetM(prhs[0]); //获取矩阵行数
	ncolsU1 = mxGetN(prhs[0]); //获取矩阵列数
	inMatU1 = mxGetPr(prhs[0]);//获取输入矩阵的指针

	//第二个输入变量的传递（xinew）
	size_t mrowsU2;    //行数
	size_t ncolsU2;    //列数
	double *inMatU2;
	mrowsU2 = mxGetM(prhs[1]); //获取矩阵行数
	ncolsU2 = mxGetN(prhs[1]); //获取矩阵列数
	inMatU2 = mxGetPr(prhs[1]);//获取输入矩阵的指针
	/* 
	这里需要注意的是Matlab中矩阵的储存是列优先的，而C语言中是行优先的，在调用矩阵元素时需要注意：
	double result;
	 将iMat中的第 i行 j列的元素值赋给result 
	result = inMat[j*mrows+i]
	*/
	//第三个输入变量的传递（网络邻接矩阵）
	size_t mrowsadj;    //行数
	size_t ncolsadj;    //列数
	double *inMatadj;  //接收输入参数的指针

	mrowsadj = mxGetM(prhs[2]); //获取矩阵行数
	ncolsadj = mxGetN(prhs[2]); //获取矩阵列数
	inMatadj = mxGetPr(prhs[2]);//获取输入矩阵的指针

	//第四个输入变量的传递（社区数）
	double numcommunity;
	numcommunity = mxGetScalar(prhs[3]);

	//第五个输入变量的传递（differindex）
	size_t mrowsdiffer;    //行数
	size_t ncolsdiffer;    //列数
	double *inMatdiffer;  //接收输入参数的指针
	mrowsdiffer = mxGetM(prhs[4]); //获取矩阵行数
	ncolsdiffer = mxGetN(prhs[4]); //获取矩阵列数
	inMatdiffer = mxGetPr(prhs[4]);//获取输入矩阵的指针
	//第六个输入变量的传递（Qg）
	double Q;
	Q = mxGetScalar(prhs[5]);

	//第七个输入变量的传递(m2)
	double m22;
	m22=mxGetScalar(prhs[6]);

	//第八个输入变量的传递(B)
	size_t mrowsB;    //行数
	size_t ncolsB;    //列数
	double *inMatB;  //接收输入参数的指针

	mrowsB = mxGetM(prhs[7]); //获取矩阵行数
	ncolsB = mxGetN(prhs[7]); //获取矩阵列数
	inMatB = mxGetPr(prhs[7]);//获取输入矩阵的指针

	size_t n=mrowsadj;
	size_t diffnum=ncolsdiffer ;
	size_t m2=int(m22);
	int c=int(numcommunity);
	 
    double *deltaQg;
    plhs[0]=mxCreateDoubleMatrix(1,diffnum,mxREAL);
	deltaQg=mxGetPr(plhs[0]);
	
	double qgTemp=Q;

	double* u0=(double*)mxMalloc(c*ncolsB*sizeof(double));
	for(size_t i=0;i<diffnum;i++){
		for(size_t k=0;k<c;k++){
			int iii=int(inMatdiffer[i*1+0]);
			int ii=iii-1;
			u0[ii*c+k]=0.0;
			for(size_t j=0;j<n;j++){
				if (j!=ii) {
					u0[ii*c+k]=u0[ii*c+k]+inMatB[j*mrowsB+ii]*inMatU1[j*mrowsU1+k];
				}
			}
//             double bbb=inMatB[ii*mrowsB+ii];
//             if(bbb==0){
//                 u0[ii*c+k]=0;
//             }
//             else{
            u0[ii*c+k]=(1/(-inMatB[ii*mrowsB+ii]))*u0[ii*c+k];
//             }
		}
	}
//     printf("%f",qcTemp);
//	double* deltaQi=(double*)mxMalloc(diffnum*sizeof(double));
	for(size_t i=0;i<diffnum;i++){
		
		int iii= int(inMatdiffer[i*1+0]);
        int ii=iii-1;
		for(int k=0;k<c;k++){
			deltaQg[i]=deltaQg[i]+((inMatU2[ii*mrowsU2+k]-inMatU1[ii*mrowsU1+k])*(inMatU2[ii*mrowsU2+k]+inMatU1[ii*mrowsU1+k]-2*u0[ii*c+k]));
		}
		deltaQg[i]=deltaQg[i]*(inMatB[ii*mrowsB+ii]/m2);
	}
	//*deltaQg=*deltaQi;
	mxFree(u0);
	u0=NULL;   
	//mxFree(deltaQi);
	//deltaQi=NULL;

}




