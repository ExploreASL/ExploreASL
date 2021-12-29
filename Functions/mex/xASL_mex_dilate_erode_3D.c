/*
 * xASL_mex_dilate_erode_3D.c
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
           double  *oima)
{
    int i=0, j=0, k=0;
    int  kerneli=0;
    int  kernelj=0;
    int  kernelk=0;
    double   kv=0.0;
    double   ov=0.0;
    double *poima;
    double *piima;
    int offset;
    double *piimakrnl;
    double *piimakrnlprep;
    double *pkrnl;
    int *offsetkernel;
    double *nonzerokernel;
    int nonzerocount;
    int *poffsetkernel;
    int  shx = (kdim[0]-1)/2;
    int  shy = (kdim[1]-1)/2;
    int  shz = (kdim[2]-1)/2;
    int shyX;
    int shzXY;
    int nxy;
    int nkdim;

    nxy = dim[0]*dim[1];
    nkdim = kdim[0]*kdim[1]*kdim[2];
    shyX  = shy*dim[0];
    shzXY = shz*dim[0]*dim[1];
    
    nonzerokernel = malloc(nkdim*sizeof(double));
    offsetkernel = malloc(nkdim*sizeof(int));
    
    /*mexPrintf("Kernel radius: %d,%d,%d\n",(int) shx,(int) shy,(int) shz);
    mexPrintf("Size: %d,%d,%d\n",(int) dim[0],(int) dim[1],(int) dim[2]);*/

    /* Initialize the input and output pointer */
    piima = iima;
    poima = oima;
    offset = 0;
    
    nonzerocount = 0;
    poffsetkernel = offsetkernel;
    pkrnl = krnl;
    /* Runs along the kernel */
    for (kernelk = -shz;kernelk <= shz; kernelk++) {
        for (kernelj = -shy;kernelj <= shy; kernelj++) {
            for (kerneli = -shx;kerneli <= shx; kerneli++) {
                /* Separate zero components of the kernel */
                if (*pkrnl>0) {
                    nonzerokernel[nonzerocount] = *pkrnl;
                    *poffsetkernel = kerneli + kernelj*dim[0] + kernelk*nxy;
                    poffsetkernel++;
                    nonzerocount++;
                }
                pkrnl++;
            }
        }
    }
    
    
    if (dilate) {
        /* Runs across the whole image */
        offset += shzXY;
        for (k=shz; k<dim[2]-shz; k++) {
            offset += shyX;
            for (j=shy; j<dim[1]-shy; j++) {
                /* Increment the image pointers to skip the gap */
                offset += shx;
                piima = iima + offset;
                poima = oima + offset;
                /* Runs along the first dimension, avoid border regions for the kernel */
                for (i=shx; i<(dim[0]-shx); i++) {
                    /* Initializes the kernel pointer */
                    pkrnl = nonzerokernel;/*krnl;*/
                    ov = *poima;

                    poffsetkernel = offsetkernel;
                    
                    for (kerneli = 0;kerneli < nonzerocount;kerneli++) {
                        piimakrnl = piima + (*poffsetkernel++);
                        kv = (*pkrnl++)*(*piimakrnl);
                        if (kv>ov)
                            ov = kv;
                    }
                                          
                    *poima = ov;
                    /* Increments the image pointers */
                    piima++;
                    poima++;
                    offset++;
                }
                offset += shx;
            }
            offset += shyX;
        }
    }
    else{
          /* Runs across the whole image */
        offset += shzXY;
        for (k=shz; k<dim[2]-shz; k++) {
            offset += shyX;
            for (j=shy; j<dim[1]-shy; j++) {
                /* Increment the image pointers to skip the gap */
                offset += shx;
                piima = iima + offset;
                poima = oima + offset;
                /* Runs along the first dimension, avoid border regions for the kernel */
                for (i=shx; i<(dim[0]-shx); i++) {
                    /* Initializes the kernel pointer */
                    pkrnl = nonzerokernel;/*krnl;*/
                    ov = *poima;

                    poffsetkernel = offsetkernel;
                    
                    for (kerneli = 0;kerneli < nonzerocount;kerneli++) {
                        piimakrnl = piima + (*poffsetkernel++);
                        kv = (*pkrnl++)*(*piimakrnl);
                        if (kv<ov)
                            ov = kv;
                    }
                                          
                    *poima = ov;
                    /* Increments the image pointers */
                    piima++;
                    poima++;
                    offset++;
                }
                offset += shx;
            }
            offset += shyX;
        }
    }
        
   free(nonzerokernel);
   free(offsetkernel);
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

    if (nrhs < 3) mexErrMsgTxt("xASL_mex_dilate_erode_3D: Not enough input arguments.");
    if (nrhs > 3) mexErrMsgTxt("xASL_mex_dilate_erode_3D: Too many input arguments.");
    if (nlhs < 1) mexErrMsgTxt("xASL_mex_dilate_erode_3D: Not enough output arguments");
    if (nlhs > 1) mexErrMsgTxt("xASL_mex_dilate_erode_3D: Too many output arguments.");

    /* Get image */

    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
    {
        mexErrMsgTxt("xASL_mex_dilate_erode_3D: ima must be numeric, real, full and double");
    }
    ndim = mxGetNumberOfDimensions(prhs[0]);
    if ((ndim < 2) || (ndim > 3))
    {
        mexErrMsgTxt("xASL_mex_dilate_erode_3D: ima must be 2- or 3-dimensional");
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

    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
    {
        mexErrMsgTxt("xASL_mex_dilate_erode_3D: kernel must be numeric, real, full and double");
    }
    krn_ndim = mxGetNumberOfDimensions(prhs[1]);
    /*if (krn_ndim != mxGetNumberOfDimensions(prhs[0]))
    {
        mexErrMsgTxt("xASL_mex_dilate_erode_3D: ima and kernel must have same dimensionality");
    }*/
    

    /* Check if dilate or erode */

    if (!mxIsChar(prhs[2]) || (mxGetM(prhs[2]) != 1))
    {
        mexErrMsgTxt("xASL_mex_dilate_erode_3D: third parameter should be a string");
    }
   
    buflen = mxGetN(prhs[2])*sizeof(mxChar)+1;
    fnc_str = (char *) mxMalloc(buflen);
    mxGetString(prhs[2],fnc_str,buflen);

    if (strcmp(fnc_str,"dilate") && strcmp(fnc_str,"erode"))
    {
        mexErrMsgTxt("xASL_mex_dilate_erode_3D: function must have value 'dilate' or 'erode'");
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
    if (krn_ndim==2) {kdim[2]=1; krn_ndim=3;} else {kdim[2]=krn_cdim[2];} 
    
    /* Execute dilation or erosion */
    do_it(iima,dim,krnl,kdim,dilate,oima);
    
    mxFree(fnc_str);
}
