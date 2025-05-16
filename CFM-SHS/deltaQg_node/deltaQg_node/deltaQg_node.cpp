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
	//nlhs�������������
	//plhs����������б�
	//nrhs�������������
	//prhs����������б�   

	// Check for proper number of input and output arguments 
	if (nrhs != 8) {
		mexErrMsgIdAndTxt("MATLAB:mxcreatecellmatrix:nrhs",
			"���������������ȷ");
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("MATLAB:mxcreatecellmatrix:nlhs",
			"���������������ȷ");
	}

	//��һ����������Ĵ��ݣ�ԭ����Ubefore��
	size_t mrowsU1;    //����
	size_t ncolsU1;    //����
	double *inMatU1;
	mrowsU1 = mxGetM(prhs[0]); //��ȡ��������
	ncolsU1 = mxGetN(prhs[0]); //��ȡ��������
	inMatU1 = mxGetPr(prhs[0]);//��ȡ��������ָ��

	//�ڶ�����������Ĵ��ݣ�xinew��
	size_t mrowsU2;    //����
	size_t ncolsU2;    //����
	double *inMatU2;
	mrowsU2 = mxGetM(prhs[1]); //��ȡ��������
	ncolsU2 = mxGetN(prhs[1]); //��ȡ��������
	inMatU2 = mxGetPr(prhs[1]);//��ȡ��������ָ��
	/* 
	������Ҫע�����Matlab�о���Ĵ����������ȵģ���C�������������ȵģ��ڵ��þ���Ԫ��ʱ��Ҫע�⣺
	double result;
	 ��iMat�еĵ� i�� j�е�Ԫ��ֵ����result 
	result = inMat[j*mrows+i]
	*/
	//��������������Ĵ��ݣ������ڽӾ���
	size_t mrowsadj;    //����
	size_t ncolsadj;    //����
	double *inMatadj;  //�������������ָ��

	mrowsadj = mxGetM(prhs[2]); //��ȡ��������
	ncolsadj = mxGetN(prhs[2]); //��ȡ��������
	inMatadj = mxGetPr(prhs[2]);//��ȡ��������ָ��

	//���ĸ���������Ĵ��ݣ���������
	double numcommunity;
	numcommunity = mxGetScalar(prhs[3]);

	//�������������Ĵ��ݣ�differindex��
	size_t mrowsdiffer;    //����
	size_t ncolsdiffer;    //����
	double *inMatdiffer;  //�������������ָ��
	mrowsdiffer = mxGetM(prhs[4]); //��ȡ��������
	ncolsdiffer = mxGetN(prhs[4]); //��ȡ��������
	inMatdiffer = mxGetPr(prhs[4]);//��ȡ��������ָ��
	//��������������Ĵ��ݣ�Qg��
	double Q;
	Q = mxGetScalar(prhs[5]);

	//���߸���������Ĵ���(m2)
	double m22;
	m22=mxGetScalar(prhs[6]);

	//�ڰ˸���������Ĵ���(B)
	size_t mrowsB;    //����
	size_t ncolsB;    //����
	double *inMatB;  //�������������ָ��

	mrowsB = mxGetM(prhs[7]); //��ȡ��������
	ncolsB = mxGetN(prhs[7]); //��ȡ��������
	inMatB = mxGetPr(prhs[7]);//��ȡ��������ָ��

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




