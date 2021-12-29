/*function imHist = xASL_im_JointHist(imA,imB,imMask,minA,maxA,minB,maxB,nBins)
* xASL_im_JointHist calculates a joint histogram of two images across a binary mask
*
* FORMAT: imHist = xASL_im_JointHist(imA,imB,imMask,minA,maxA,minB,maxB,nBins)
*
* INPUT:
*   imA    - First input image (REQUIRED).
*   imB    - Second input image, needs to have the same dimensions (REQUIRED).
*   imMask - Binary mask to calculate the joint histogram (REQUIRED).
*   minA   - Minimal value to be counted for image A (REQUIRED).
*   maxA   - Maximal value to be counted for image A (REQUIRED).
*   minB   - Minimal value to be counted for image B (REQUIRED).
*   minB   - Maximal value to be counted for image B (REQUIRED).
*   nBins  - Number of bins (REQUIRED).
*
* OUTPUT:
*   imHist - Resulting joint histogram of size nBins x nBins.
*
* -----------------------------------------------------------------------------------------------------------------------------------------------------
* DESCRIPTION: It calculates a joint histogram of two images of any dimensions over a mask of the same size.
*              Values outside of the bins are counted to the first/last bin.
* EXAMPLE: 
*     imHist = xASL_im_JointHist(imA,imB,imMask,0,200,-10,100,50)
* -----------------------------------------------------------------------------------------------------------------------------------------------------
* __________________________________
* Copyright © 2015-2019 ExploreASL
*
* 2019-06-26 JP*/

#define myCreateMatrix(a,b) 

#include <string.h>
#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	FR_IN	  prhs[0]
#define	FT_IN	  prhs[1]
#define FM_IN     prhs[2]
#define MIN1_IN   prhs[3]
#define MAX1_IN   prhs[4]
#define MIN2_IN   prhs[5]
#define MAX2_IN   prhs[6]
#define BINS_IN   prhs[7]

/* Output Arguments */

#define	JH_OUT	plhs[0]

#define NAME "xASL_mex_JointHist"

#ifndef MAX
#define MAX(A,B) ((A>B) ? (A) : (B))
#endif
#ifndef MIN
#define MIN(A,B) ((A<B) ? (A) : (B))
#endif

/*xASL_mex_JointHist(im1,im2,imMask,min1,max1,min2,max2,bins)*/
void mexFunction(
		int nlhs,       mxArray *plhs[],
		int nrhs, const mxArray *prhs[]
		) {
	int numel,bins,i,xf,yf;
	double y,x,min1,max1,min2,max2,*fr,*ft,*JH,*fm;
	
	if ((nrhs != 8)) {
		mexErrMsgTxt(NAME " requires eight input arguments.");
	} else if ((nlhs != 1)) {
		mexErrMsgTxt(NAME " requires one output argument.");
	}
 
	for (i=0;i<nrhs;i++) {
		if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) ||
				mxIsSparse(prhs[i])  || !mxIsDouble(prhs[i]))
			mexErrMsgTxt(NAME " needs all input arguments to be real.");
	}
	
	numel = mxGetNumberOfElements(FR_IN);
	fr    = mxGetPr(FR_IN);
	ft    = mxGetPr(FT_IN);
	fm    = mxGetPr(FM_IN);
	min1  = (double)(*(mxGetPr(MIN1_IN))) ;
	max1  = (double)(*(mxGetPr(MAX1_IN))) ;
	min2  = (double)(*(mxGetPr(MIN2_IN))) ;
	max2  = (double)(*(mxGetPr(MAX2_IN))) ;
	bins  = (int)(*(mxGetPr(BINS_IN)));
	
	JH_OUT = mxCreateDoubleMatrix(bins,bins,mxREAL); 
	JH     = mxGetPr(JH_OUT) ;
	
	for (i=0;i<bins*bins;i++) {
		JH[i] = 0;
	}
	
	for (i=0;i<numel;i++) {
		if (fm[i]) {
			y = (fr[i]-min1)/(max1-min1)*bins;
			x = (ft[i]-min2)/(max2-min2)*bins;
			
			yf = (int)(round(y));
			xf = (int)(round(x));
			
			yf = MAX(yf,0);
			xf = MAX(xf,0);
			
			yf = MIN(yf,bins-1);
			xf = MIN(xf,bins-1);
			
			JH[yf + xf*bins] += 1;
			
		}
	}
	
}
