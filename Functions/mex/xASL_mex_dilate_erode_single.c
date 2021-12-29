/*
 * xASL_mex_dilate_erode_single.c
 * Jan Petr
 */

#include <string.h>
#include "mex.h"

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

/* Function doing the job. */

void do_it(double  *iima,
           mwSize  dim[3],
           double  *krnl,
           mwSize  kdim[3],
           int     dilate,
           double  *oima,
	   int     dimdir)
{
    int  i=0, j=0, k=0;
    int  kerneli=0;
    double   kv=0.0;
    double   ov=0.0;
    double *poima;
    double *piima;
    double *piimakrnl;
    double *pkrnl;
    int  sh = (kdim[0]-1)/2;
    int shx;
    int shxy;
    int nxy;

    nxy = dim[0]*dim[1];
    shx = sh*dim[0];
    shxy = sh*dim[0]*dim[1];
    
    /*mexPrintf("Kernel radius: %d, direction: %d\n",(int) sh,dimdir);*/

    /* Initialize the input and output pointer */
    piima = iima;
    poima = oima;
    
    switch (dimdir)
    {
        case 1:
            /* Runs across the whole image */
            for (k=0; k<dim[2]; k++) {
                for (j=0; j<dim[1]; j++) {
                    /* Increment the image pointers to skip the gap */
                    piima += sh;
                    poima += sh;
                    /* Runs along the first dimension, avoid border regions for the kernel */
                    for (i=sh; i<(dim[0]-sh); i++) {
                        /* Initializes the kernel pointer */
                        pkrnl = krnl;
                        ov = *poima;
                        piimakrnl = piima-sh;
                        /* Runs along the kernel */
                        for (kerneli = -sh;kerneli <= sh; kerneli++) {
                            kv = (*pkrnl++)*(*piimakrnl++);
                            if (dilate)
                                ov = MAX(ov,kv);
                            else
                                ov = MIN(ov,kv);
                        }
                        *poima = ov;
                        /* Increments the image pointers */
                        piima++;
                        poima++;
                    }
                    piima += sh;
                    poima+=sh;
                }
            }
            break;
        case 2:
            for (k=0; k<dim[2]; k++) {
                /* Increment the image pointers to skip the gap */
                piima += shx;
                poima += shx;
                for (j=sh; j<(dim[1]-sh); j++) {
                    for (i=0; i<dim[0]; i++) {
                        pkrnl = krnl;
                        ov = *poima;
                        piimakrnl = piima-shx;
                        for (kerneli = -sh;kerneli <= sh; kerneli++) {
                            kv = (*pkrnl++)*(*piimakrnl);
                            piimakrnl += dim[0];
                            if (dilate)
                                ov = MAX(ov,kv);
                            else
                                ov = MIN(ov,kv);
                        }
                        *poima = ov;
                        piima++;
                        poima++;
                    }
                }
                piima += shx;
                poima += shx;
            }
            break;
        case 3:
            /* Increment the image pointers to skip the gap */
            piima += shxy;
            poima += shxy;
            for (k=sh; k<(dim[2]-sh); k++) {
                for (j=0; j<dim[1]; j++) {
                    for (i=0; i<dim[0]; i++) {
                        pkrnl = krnl;
                        ov = *poima;
                        piimakrnl = piima-shxy;
                        for (kerneli = -sh;kerneli <= sh; kerneli++) {
                            kv = (*pkrnl++)*(*piimakrnl);
                            piimakrnl += nxy;
                            if (dilate)
                                ov = MAX(ov,kv);
                            else
                                ov = MIN(ov,kv);
                        }
                        *poima = ov;
                        piima++;
                        poima++;
                    }
                }
            }
            break;
    }
}


/* Gateway function */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char          *fnc_str = NULL;
    mwSize        ndim=0, krn_ndim=0;
    int           n, i;
    int           dilate = 0;
    mwSize        buflen = 0;
    const mwSize  *cdim = NULL, *krn_cdim = NULL;
    mwSize        dim[3], kdim[3];
    double        *iima = NULL;
    double        *oima = NULL;
    double        *krnl = NULL;
    double        *krnly = NULL;
    double        *krnlz = NULL;

    if (nrhs < 5) mexErrMsgTxt("xASL_mex_dilate_erode_single: Not enough input arguments.");
    if (nrhs > 5) mexErrMsgTxt("xASL_mex_dilate_erode_single: Too many input arguments.");
    if (nlhs < 1) mexErrMsgTxt("xASL_mex_dilate_erode_single: Not enough output arguments");
    if (nlhs > 1) mexErrMsgTxt("xASL_mex_dilate_erode_single: Too many output arguments.");

    /* Get image */

    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
    {
        mexErrMsgTxt("xASL_mex_dilate_erode_single: ima must be numeric, real, full and double");
    }
    ndim = mxGetNumberOfDimensions(prhs[0]);
    if ((ndim < 2) || (ndim > 3))
    {
        mexErrMsgTxt("xASL_mex_dilate_erode_single: ima must be 2- or 3-dimensional");
    }
    cdim = mxGetDimensions(prhs[0]);
    iima = mxGetPr(prhs[0]);

    /* specifies along which dimension to operate */
    /*dimdir=(int)(*(mxGetPr(prhs[3])));*/

    /* Fix dimensions to allow for 2D and 3D data */

    dim[0] = cdim[0]; dim[1] = cdim[1];
    if (ndim==2) {dim[2]=1; ndim=3;} else {dim[2]=cdim[2];} 

    /* Calculate the number of voxels */
    for (i=0, n=1; i<ndim; i++)
    {
        n *= dim[i];
    }

    /* Get kernel */

    for (i=1;i<4;i++) {
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) || !mxIsDouble(prhs[i]))
        {
            mexErrMsgTxt("xASL_mex_dilate_erode_single: kernel must be numeric, real, full and double");
        }
        krn_ndim = mxGetNumberOfDimensions(prhs[i]);
        if (krn_ndim > 2)
        {
            mexErrMsgTxt("xASL_mex_dilate_erode_single: kernel must be a vector Nx1");
        }
    }

    /* Check if dilate or erode */

    if (!mxIsChar(prhs[4]) || (mxGetM(prhs[4]) != 1))
    {
        mexErrMsgTxt("xASL_mex_dilate_erode_single: third parameter should be a string");
    }
   
    buflen = mxGetN(prhs[4])*sizeof(mxChar)+1;
    fnc_str = (char *) mxMalloc(buflen);
    mxGetString(prhs[4],fnc_str,buflen);

    if (strcmp(fnc_str,"dilate") && strcmp(fnc_str,"erode"))
    {
        mexErrMsgTxt("xASL_mex_dilate_erode_single: function must have value 'dilate' or 'erode'");
    }
    if (!strcmp(fnc_str,"dilate")) {dilate = 1;}
    else {dilate = 0;}
   
    /* Allocate and initialise output image */

    plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
                                   mxGetDimensions(prhs[0]),mxDOUBLE_CLASS,mxREAL);
    oima = mxGetPr(plhs[0]);
    memcpy(oima,iima,n*sizeof(double));

    /* Load the kernel */
    krn_cdim = mxGetDimensions(prhs[1]);
    krnl = mxGetPr(prhs[1]);
    kdim[0]=krn_cdim[0]; kdim[1]=krn_cdim[1];
    
    /* Skip too short dimensions */
    if (kdim[0]>1) {
    /* Execute dilation or erosion */
    do_it(iima,dim,krnl,kdim,dilate,oima,1);
    
    /* Copy results to source and repeat for the other two dimensions */
    memcpy(iima,oima,n*sizeof(double));
    }
    
    krn_cdim = mxGetDimensions(prhs[2]);
    krnl = mxGetPr(prhs[2]);
    kdim[0]=krn_cdim[0]; kdim[1]=krn_cdim[1];
    
    if (kdim[0]>1) {
    do_it(iima,dim,krnl,kdim,dilate,oima,2);
    memcpy(iima,oima,n*sizeof(double));
    }
    
    krn_cdim = mxGetDimensions(prhs[3]);
    krnl = mxGetPr(prhs[3]);
    kdim[0]=krn_cdim[0]; kdim[1]=krn_cdim[1];
    
    if (kdim[0]>1) {
    do_it(iima,dim,krnl,kdim,dilate,oima,3);
    }
    mxFree(fnc_str);
}
