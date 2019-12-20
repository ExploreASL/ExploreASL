/*
 * xASL_mex_conv3DsepGauss.c
 * 3D separable convolution with a Gaussian kernel
 * Float is on the input 
 * Returned is the convoluted image
 * On input is the sigma in 3D
 *
 * [imConv] = xASL_mex_conv3Dsep(im,[X,Y,Z])
 * im - 3D image, double
 * [X,Y,Z] - sigma in voxels of the separable 3D Gaussian kernel in all three dimensions
 *         - put zero and this dimension will not be convoluted
 * imConv - convolved image
 * 
 * Implemented by Jan Petr
 */

#include <string.h>
#include <math.h>
#include "mex.h"

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
		/* The central point has them all */
		for (i=0;i<=nkx;i++) kernelNorm[0] += kernelX[i+nkx];
		for (i=-nkx;i<=0;i++) kernelNormEnd[0] += kernelX[i+nkx];
		
		/* Enlarge the sum */
		for (i=1;i<=nkx;i++)  kernelNorm[i] = kernelNorm[i-1] + kernelX[nkx-i];
		for (i=1;i<=nkx;i++)  kernelNormEnd[i] = kernelNormEnd[i-1] + kernelX[nkx+i];
		
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
		/* The central point has them all */
		for (i=0;i<=nky;i++) kernelNorm[0] += kernelY[i+nky];
		for (i=-nky;i<=0;i++) kernelNormEnd[0] += kernelY[i+nky];
		/* Enlarge the sum */
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
    mwSize        dim[3];
	int           window[3];
	double        *kernel[3];	
	double        kernelSum;
    double        *iima = NULL;
    double        *oima = NULL;
	double        *wima = NULL;
	double        *sigma = NULL;

    if (nrhs < 2) mexErrMsgTxt("xASL_mex_conv3DsepGauss: Not enough input arguments.");
    if (nrhs > 2) mexErrMsgTxt("xASL_mex_conv3DsepGauss: Too many input arguments.");
    if (nlhs < 1) mexErrMsgTxt("xASL_mex_conv3DsepGauss: Not enough output arguments");
    if (nlhs > 1) mexErrMsgTxt("xASL_mex_conv3DsepGauss: Too many output arguments.");

    /* Get image */

    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
    {
        mexErrMsgTxt("xASL_mex_conv3DsepGauss: ima must be numeric, real, full and double");
    }
	
	if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
    {
        mexErrMsgTxt("xASL_mex_conv3DsepGauss: kernel dimensions must be numeric, real, full and double");
    }
	
    ndim = mxGetNumberOfDimensions(prhs[0]);
    if ((ndim < 2) || (ndim > 3))
    {
        mexErrMsgTxt("xASL_mex_conv3DsepGauss: ima must be 2- or 3-dimensional");
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

	cdim = mxGetDimensions(prhs[1]);
	sigma = mxGetPr(prhs[1]);
	
	window[0] = (mwSize)ceil(3*sigma[0]);
	if ((sigma[0]>0) & ((2*3*sigma[0]) > dim[0]))
	{
		mexErrMsgTxt("xASL_mex_conv3DsepGauss: Kernel dimension 1 too large.");
	}
	
	if ( (cdim[0]>1) | (cdim[1]>1) )
	{
		window[1] = (mwSize)ceil(3*sigma[1]);
		if ((sigma[1]>0) & ((2*3*sigma[1]) > dim[1]))
		{
			mexErrMsgTxt("xASL_mex_conv3DsepGauss: Kernel dimension 2 too large.");
		}
	}
	else
	{
		window[1] = 0;
	}
	
	if ( (cdim[0]>2) | (cdim[1]>2) )
	{
		window[2] = (mwSize)ceil(3*sigma[2]);
		if ((sigma[2]>0) & ((2*3*sigma[2]) > dim[2]))
		{
			mexErrMsgTxt("xASL_mex_conv3DsepGauss: Kernel dimension 3 too large.");
		}
	}
	else
	{
		window[2] = 0;
	}
		
	/* Allocate the three kernels */
	
	for (i=0;i<3;i++)
	{
		if (window[i]>0)
		{
			kernel[i] = (double *) mxMalloc(sizeof(double)*(2*window[i]+1));
			
			kernelSum = 0;
			for (k=-(window[i]);k<=window[i];k++)
			{
				kernel[i][k+window[i]] = exp(((double)-k)*((double)k)/(2.0*sigma[i]*sigma[i]));
				kernelSum += kernel[i][k+window[i]];
			}
			if (kernelSum)
			{
				for (k=-window[i];k<=window[i];k++)
				{
					kernel[i][k+window[i]] = kernel[i][k+window[i]]/kernelSum;
					
					
				}
			}
		}
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
    do_it(wima,dim,oima,window,kernel[0],kernel[1],kernel[2]);
    
	for (i=0;i<3;i++)
	{
		if (window[i]>0)
			mxFree(kernel[i]);
	}
	
	mxFree(wima);

}
