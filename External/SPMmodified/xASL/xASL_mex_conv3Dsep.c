/*
 * xASL_mex_conv3Dsep.c
 * 3D separable convolution with a supplied kernel
 * Float is on the input 
 * Returned is the convoluted image
 *
 * [imConv] = xASL_mex_conv3Dsep(im,kX,[kY,kZ])
 * im - 3D image, double
 * kX, kY, kZ - are the 1D kernels for the convolution
 * 	- kY and kZ are optional
 *      - if 0 or [] is supplied, then convolution in this dimension will be skipped
 * imConv - convolved image
 * 
 * Implemented by Jan Petr
 */

#include <string.h>
#include <math.h>
#include "mex.h"

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif

/* Function doing the job. */

void do_it(double  *iima,
           mwSize  dim[3],
           double  *oima,
		   int     *window,
		   double  *kernelX,
           double  *kernelY,
           double  *kernelZ)
{
    int i,j,k,z;
    int n,nx,nxy,nyz,nz,ny;
	int nkx,nky,nkz,nkMar;
    double *piima;
    double *poima;
	double *kernelNorm;
	double *kernelNormEnd;
	double *pkernel;
	double *piimaloc;
	double sum;
    
    n = dim[0]*dim[1]*dim[2];
    nx = dim[0];
    nxy = dim[0]*dim[1];
	nyz = dim[1]*dim[2];
	nz = dim[2];
	ny = dim[1];
	nkx = window[0];
	nky = window[1];
	nkz = window[2];

	if (nkx)
	{
		/******* Dimension X *********/
		kernelNorm = (double *) mxMalloc(sizeof(double)*(nkx+1));
		kernelNormEnd = (double *) mxMalloc(sizeof(double)*(nkx+1));
		kernelNorm[0] = 0;
		kernelNormEnd[0] = 0;
		/* Calculates the norm to divide the weighted sum */
		/* For border points, we get the center of the kernel and one half of it */
		for (i=0;i<=nkx;i++) kernelNorm[0] += kernelX[i+nkx];
		for (i=-nkx;i<=0;i++) kernelNormEnd[0] += kernelX[i+nkx];
		
		/* Once we move from the border, then the part of the kernel over which we convolve gets bigger */
		for (i=1;i<=nkx;i++)  kernelNorm[i] = kernelNorm[i-1] + kernelX[nkx-i];
		for (i=1;i<=nkx;i++)  kernelNormEnd[i] = kernelNormEnd[i-1] + kernelX[nkx+i];
		
		/* Sometimes, the dimension is smaller than the kernel - then we cannot filter over the whole
		 * kernel, not even in the middle of the dimension */
		nkMar = MIN(nkx,nx-nkx);
		
		/* In case the kernel is larger than the dimension, then we need to start subtracting - only from the starting
		 * one, the endnorm is not going to be used*/
		if (nx < (2*nkx-1))
		{
			/* For all those missing, we can calculate the norm with cutting the filter at both ends */
			for (i=nkMar;i<nx-nkMar;i++)
			{
				for (j=nkx-i+nx;j<=(2*nkx);j++) kernelNorm[i] -= kernelX[j];
			}
		}
		/* Go across Y*Z */
		poima = oima;
		piima = iima;
		for (j=0;j<nyz;j++)
		{
			/* Intro margin */
			for (i=0;i<nkMar;i++)
			{
				sum = 0;
				pkernel = kernelX-i+nkx;
				piimaloc = piima - i;
				/* Go across the kernel */
				for (k=-i;k<=nkx;k++)
					sum += (*piimaloc++)*(*pkernel++);
				
				*poima++ = sum/kernelNorm[i];
				piima++;
			}
			
			/* In case the kernel does not fit entirely, we have to do an extra filtering of the center */
			if (nx < (2*nkx-1))
			{
				for (i=nkMar;i<nx-nkMar;i++)
				{
					sum = 0;
					pkernel = kernelX - i + nkx;
					piimaloc = piima - i;
					// Go across the kernel
					for (k=-i;k<nx-i;k++)
					{
						sum += (*piimaloc++)*(*pkernel++);
					}
					
					*poima++ = sum/kernelNorm[i];
					piima ++;
				}
			}
			else
			{
				/* Middle of the first dimension */
				for (i=nkx;i<(nx-nkx);i++)
				{
					sum = 0;
					pkernel = kernelX;
					piimaloc = piima - nkx;
					/* Go across the kernel */
					for (k=-nkx;k<=nkx;k++)
						sum += (*piimaloc++)*(*pkernel++);
					
					*poima++ = sum;
					piima++;
				}
			}
			/* End margin */
			for (i=nx-nkMar;i<nx;i++)
			{
				sum = 0;
				/* Go across the kernel */
				pkernel = kernelX;
				piimaloc = piima-nkx;
				for (k=-nkx;k<(nx-i);k++)
					sum += (*piimaloc++)*(*pkernel++);
				
				*poima++ = sum/kernelNormEnd[nx-i-1];
				piima++;
			}
			
		}
		mxFree(kernelNorm);
		mxFree(kernelNormEnd);
		
		/* COPY OUTPUT TO INPUT */
		piima = iima;
		poima = oima;
		for (i=0;i<n;i++)
			*piima++ = *poima++;
	}
	
	if (nky)
	{
		/******* Dimension Y *********/
		kernelNorm = (double *) mxMalloc(sizeof(double)*(nky+1));
		kernelNormEnd = (double *) mxMalloc(sizeof(double)*(nky+1));
		kernelNorm[0] = 0;
		kernelNormEnd[0] = 0;
		/* Calculates the norm to divide the weighted sum */
		/* For border points, we get the center of the kernel and one half of it */
		for (i=0;i<=nky;i++) kernelNorm[0] += kernelY[i+nky];
		for (i=-nky;i<=0;i++) kernelNormEnd[0] += kernelY[i+nky];
		
		/* Once we move from the border, then the part of the kernel over which we convolve gets bigger */
		for (i=1;i<=nky;i++)  kernelNorm[i] = kernelNorm[i-1] + kernelY[nky-i];
		for (i=1;i<=nky;i++)  kernelNormEnd[i] = kernelNormEnd[i-1] + kernelY[nky+i];
		
		/* Sometimes, the dimension is smaller than the kernel - then we cannot filter over the whole
		 * kernel, not even in the middle of the dimension */
		nkMar = MIN(nky,ny-nky);
		
		/* In case the kernel is larger than the dimension, then we need to start subtracting - only from the starting
		 * one, the endnorm is not going to be used*/
		if (ny < (2*nky-1))
		{
			/* For all those missing, we can calculate the norm with cutting the filter at both ends */
			for (i=nkMar;i<ny-nkMar;i++)
			{
				for (j=nky-i+ny;j<=(2*nky);j++) kernelNorm[i] -= kernelY[j];
			}
		}
		
		/* Go across X*Z */
		for (z=0;z<nz;z++)
		{
			for (j=0;j<nx;j++)
			{
				poima = oima + j + z*nxy;
				piima = iima + j + z*nxy;
				/* Intro margin */
				for (i=0;i<nkMar;i++)
				{
					sum = 0;
					pkernel = kernelY-i+nky;
					piimaloc = piima - i*nx;
					/* Go across the kernel */
					for (k=-i;k<=nky;k++)
					{
						sum += (*piimaloc)*(*pkernel++);
						piimaloc += nx;
					}
					
					*poima = sum/kernelNorm[i];
					poima += nx;
					piima += nx;
				}
				
				/* In case the kernel does not fit entirely, we have to do an extra filtering of the center */
				if (ny < (2*nky-1))
				{
					for (i=nkMar;i<ny-nkMar;i++)
					{
						sum = 0;
						pkernel = kernelY - i + nky;
						piimaloc = piima - i*nx;
						// Go across the kernel
						for (k=-i;k<ny-i;k++)
						{
							sum += (*piimaloc)*(*pkernel++);
							piimaloc += nx;
						}
						
						*poima = sum/kernelNorm[i];
						poima += nx;
						piima += nx;
					}
				}
				else
				{
					/* Middle of the first dimension */
					for (i=nky;i<(ny-nky);i++)
					{
						sum = 0;
						pkernel = kernelY;
						piimaloc = piima - nx*nky;
						/* Go across the kernel */
						for (k=-nky;k<=nky;k++)
						{
							sum += (*piimaloc)*(*pkernel++);
							piimaloc += nx;
						}
						
						*poima = sum;
						poima += nx;
						piima += nx;
					}
				}
				
				/* End margin */
				for (i=ny-nkMar;i<ny;i++)
				{
					sum = 0;
					/* Go across the kernel */
					pkernel = kernelY;
					piimaloc = piima - nx*nky;
					for (k=-nky;k<(ny-i);k++)
					{
						sum += (*piimaloc)*(*pkernel++);
						piimaloc += nx;
					}
					
					*poima = sum/kernelNormEnd[ny-i-1];
					poima += nx;
					piima += nx;
				}
			}
		}
		mxFree(kernelNorm);
		mxFree(kernelNormEnd);
		
		/* COPY OUTPUT TO INPUT */
		piima = iima;
		poima = oima;
		for (i=0;i<n;i++)
			*piima++ = *poima++;
	}
	
	if (nkz)
	{
		/******* Dimension Z *********/
		kernelNorm = (double *) mxMalloc(sizeof(double)*(nkz+1));
		kernelNormEnd = (double *) mxMalloc(sizeof(double)*(nkz+1));
		kernelNorm[0] = 0;
		kernelNormEnd[0] = 0;
		/* Calculates the norm to divide the weighted sum */
		/* For border points, we get the center of the kernel and one half of it */
		for (i=0;i<=nkz;i++) kernelNorm[0] += kernelZ[i+nkz];
		for (i=-nkz;i<=0;i++) kernelNormEnd[0] += kernelZ[i+nkz];
		
		/* Once we move from the border, then the part of the kernel over which we convolve gets bigger */
		for (i=1;i<=nkz;i++)  kernelNorm[i] = kernelNorm[i-1] + kernelZ[nkz-i];
		for (i=1;i<=nkz;i++)  kernelNormEnd[i] = kernelNormEnd[i-1] + kernelZ[nkz+i];
		
		/* Sometimes, the dimension is smaller than the kernel - then we cannot filter over the whole
		 * kernel, not even in the middle of the dimension */
		nkMar = MIN(nkz,nz-nkz);
		
		/* In case the kernel is larger than the dimension, then we need to start subtracting - only from the starting
		   one, the endnorm is not going to be used*/
		if (nz < (2*nkz-1))
		{
			/* For all those missing, we can calculate the norm with cutting the filter at both ends */
			for (i=nkMar;i<nz-nkMar;i++)
			{
				for (j=nkz-i+nz;j<=(2*nkz);j++) kernelNorm[i] -= kernelZ[j];
			}
		}
		
		/* Go across X*Y */
		
		for (j=0;j<nxy;j++)
		{
			poima = oima + j;
			piima = iima + j;
			// Intro margin 
			for (i=0;i<nkMar;i++)
			{
				sum = 0;
				pkernel = kernelZ - i + nkz;
				piimaloc = piima - i*nxy;
				// Go across the kernel 
				for (k=-i;k<nkz;k++)
				{
					sum += (*piimaloc)*(*pkernel++);
					piimaloc += nxy;
				}
				
				*poima = sum/kernelNorm[i];
				poima += nxy;
				piima += nxy;
			}
			/* In case the kernel does not fit entirely, we have to do an extra filtering of the center */
			if (nz < (2*nkz-1))
			{
				for (i=nkMar;i<nz-nkMar;i++)
				{
					sum = 0;
					pkernel = kernelZ - i + nkz;
					piimaloc = piima - i*nxy;
					// Go across the kernel
					for (k=-i;k<nz-i;k++)
					{
						sum += (*piimaloc)*(*pkernel++);
						piimaloc += nxy;
					}
					
					*poima = sum/kernelNorm[i];
					poima += nxy;
					piima += nxy;
				}
			}
			else
			{
				// Middle of the first dimension
				for (i=nkz;i<(nz-nkz);i++)
				{
					sum = 0;
					pkernel = kernelZ;
					piimaloc = piima - nxy*nkz;
					// Go across the kernel
					for (k=-nkz;k<nkz;k++)
					{
						sum += (*piimaloc)*(*pkernel++);
						piimaloc += nxy;
					}
					
					*poima = sum;
					poima += nxy;
					piima += nxy;
				}
			}
			
			// End margin
			for (i=nz-nkMar;i<nz;i++)
			{
				sum = 0;
				// Go across the kernel
				pkernel = kernelZ;
				piimaloc = piima - nxy*nkz;
				for (k=-nkz;k<(nz-i);k++)
				{
					sum += (*piimaloc)*(*pkernel++);
					piimaloc += nxy;
				}
				
				*poima = sum/kernelNormEnd[nz-i-1];
				poima += nxy;
				piima += nxy;
			}
		}
	
		mxFree(kernelNorm);
		mxFree(kernelNormEnd);
	}
}

/* Gateway function */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize        ndim=0;
    int           i;
	int           k;
    const mwSize  *cdim = NULL;
	int           kdim = 0;
    mwSize        dim[3];
	int           window[3];
	double        *kernelX = NULL;
	double        *kernelY = NULL;
	double        *kernelZ = NULL;
    double        *iima = NULL;
	double        *wima = NULL;
    double        *oima = NULL;

    if (nrhs < 2) mexErrMsgTxt("xASL_mex_conv3Dsep: Not enough input arguments.");
    if (nrhs > 4) mexErrMsgTxt("xASL_mex_conv3Dsep: Too many input arguments.");
    if (nlhs < 1) mexErrMsgTxt("xASL_mex_conv3Dsep: Not enough output arguments");
    if (nlhs > 1) mexErrMsgTxt("xASL_mex_conv3Dsep: Too many output arguments.");

    /* Get image */

    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
    {
        mexErrMsgTxt("xASL_mex_conv3Dsep: ima must be numeric, real, full and double");
    }
	
    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
    {
        mexErrMsgTxt("xASL_mex_conv3Dsep: kernel must be numeric, real, full and double");
    }
	
    ndim = mxGetNumberOfDimensions(prhs[0]);
    if ((ndim < 2) || (ndim > 3))
    {
        mexErrMsgTxt("xASL_mex_conv3Dsep: ima must be 2- or 3-dimensional");
    }
    cdim = mxGetDimensions(prhs[0]);
    iima = mxGetPr(prhs[0]);
   
    /* Fix dimensions to allow for 2D and 3D data */

    dim[0] = cdim[0]; dim[1] = cdim[1];
    if (ndim==2) {dim[2]=1; ndim=3;} else {dim[2]=cdim[2];} 
  
    /* Allocate and initialise output images */

    plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
                                   mxGetDimensions(prhs[0]),mxDOUBLE_CLASS,mxREAL);
    oima = mxGetPr(plhs[0]);
    
	/* First kernel should always be there */
	
	/* If the array is empty than set to skip this dimension in convolution */
	if (mxIsEmpty(prhs[1]))
	{
		window[0] = 0;
	}
	else
	{
		/* Get the dimension of the array - as maximum of the 2 first dimensions as we expect a 1D vector but can be Nx1 or 1xN */
		cdim = mxGetDimensions(prhs[1]);
		kdim = MAX(cdim[0],cdim[1]);
		
		/* Just a scalar still does not require convolution */
		if (kdim == 1)
		{
			window[0] = 0;
		}
		else
		{
			window[0] = (kdim-1)/2;
			kernelX = mxGetPr(prhs[1]);
		}
	}
	
	/* Same thing for kernelY, but it also checks the existence of the field  */
	if ( (nrhs<3) | mxIsEmpty(prhs[2]))
	{
		window[1] = 0;
	}
	else
	{
		/* Get the dimension of the array - as maximum of the 2 first dimensions as we expect a 1D vector but can be Nx1 or 1xN */
		cdim = mxGetDimensions(prhs[2]);
		kdim = MAX(cdim[0],cdim[1]);
		
		/* Just a scalar still does not require convolution */
		if (kdim == 1)
		{
			window[1] = 0;
		}
		else
		{
			window[1] = (kdim-1)/2;
			kernelY = mxGetPr(prhs[2]);
		}
	}
	
	/* Same thing for kernelZ, but it also checks the existence of the field  */
	if ( (nrhs<4) | mxIsEmpty(prhs[3]))
	{
		window[2] = 0;
	}
	else
	{
		/* Get the dimension of the array - as maximum of the 2 first dimensions as we expect a 1D vector but can be Nx1 or 1xN */
		cdim = mxGetDimensions(prhs[3]);
		kdim = MAX(cdim[0],cdim[1]);
		
		/* Just a scalar still does not require convolution */
		if (kdim == 1)
		{
			window[2] = 0;
		}
		else
		{
			window[2] = (kdim-1)/2;
			kernelZ = mxGetPr(prhs[3]);
		}
	}
	
	if ( (window[0]>0) & ( (window[0] + 2) > dim[0]))
	{
		mexErrMsgTxt("xASL_mex_conv3Dsep: Kernel dimension 1 too large.");
	}
	
	if ( (window[1]>0) & ( (window[1] + 2) > dim[1]))
	{
		mexErrMsgTxt("xASL_mex_conv3Dsep: Kernel dimension 2 too large.");
	}
	
	if ( (window[2]>0) & ( (window[2] + 2) > dim[2]))
	{
		mexErrMsgTxt("xASL_mex_conv3Dsep: Kernel dimension 3 too large.");
	}

	
	/* The algorithm needs and input and output image, and the input image is changed as well */
	/* So to avoid modify the input of this function, we need to copy the input matrix to a working */
	/* matrix and supply this to the core function */
	wima = (double *) mxMalloc(sizeof(double)*(dim[0]*dim[1]*dim[2]));
	for (i=0;i<(dim[0]*dim[1]*dim[2]);i++)
	{
		wima[i] = iima[i];
	}
	
	/* Execute the convolution */
    do_it(wima,dim,oima,window,kernelX,kernelY,kernelZ);    
	
	mxFree(wima);
}
