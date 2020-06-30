/*
 * Robert Dahnke
 * $Id: cat_vbdist.c 1523 2019-11-21 23:12:24Z gaser $ 
 *
 */

/* voxelbased euclidean distance calculation
 * ________________________________________________________________________
 * Calculates the euclidean distance without PVE to an object in P with a 
 * boundary of 0.5.
 * 
 *  [D,I,L] = vbdist(P[,R])
 *  
 *  P (single)  input image with zero for non elements and uint8 values for labels
 *  R (logical) range for distance calculation
 *  D (single)  distance image
 *  L (uint8)   label map
 *  I (uint32)  index of nearest point
 * ________________________________________________________________________
 * Robert Dahnke 2010_01
 * Center of Neuroimaging 
 * University Jena
 */

#include "mex.h"   
#include "math.h"
#include "float.h"
#include <stdio.h>
#include <stdlib.h>

/* estimate minimum of A and its index in A */
void pmin(float A[], int sA, float *minimum, int *index)
{
  int i; 
  *minimum = FLT_MAX; *index = 0; /* printf("%d ",sizeof(A)/8); */
  for(i=0; i<sA; i++) {
    if ((A[i]>0) && (*minimum>A[i]))
    { 
      *minimum = A[i]; 
      *index   = i;
    }
  }
}


/* estimate x,y,z position of index i in an array size sx,sxy=sx*sy... */
void ind2sub(int i, int *x, int *y, int *z, int sxy, int sy) {
  *z = (int)floor( (double)i / (double)sxy ) +1; 
   i = i % (sxy);
  *y = (int)floor( (double)i / (double)sy ) +1;        
  *x = i % sy + 1;
}

/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  mwSize *sL;
  int     dL;
  int     nL;
  int     x;
  int     y;
  int     xy;

  mwSize sS[2] = {1,3};
  mxArray *SS;
  double  *S; 

  float s1,s2,s3;
  float   s12;
  float   s13;
  float   s23;
  float   s123;

  int   NI[14];  
  float ND[14];
  int   sN;    
  float       *DN;
  float       DNm;
  int i, n, ni, DNi;

  float         *D;
  unsigned int  *I;
  unsigned char *L;
  
  float         *V;
  bool          *R; 
  bool e255; 
  
  int u,v,w,nu,nv,nw; 

  if (nrhs<1)                                       mexErrMsgTxt("ERROR:cat_vbdist: not enough input elements\n");
  if (nlhs>3)                                       mexErrMsgTxt("ERROR:cat_vbdist: too many output elements.\n");
  if (mxIsSingle(prhs[0])==0)                       mexErrMsgTxt("ERROR:cat_vbdist: first  input must be an 3d single matrix\n");
  if (nrhs==2 && mxIsLogical(prhs[1])==0)           mexErrMsgTxt("ERROR:cat_vbdist: second input must be an 3d logical matrix\n");
  if (nrhs==3 && mxIsDouble(prhs[2])==0)            mexErrMsgTxt("ERROR:cat_vbdist: third input must be an double matrix\n");
  if (nrhs==3 && mxGetNumberOfElements(prhs[2])!=3) mexErrMsgTxt("ERROR:cat_vbdist: third input must have 3 Elements"); 
  
  /* main information about input data (size, dimensions, ...) */
  sL = mxGetDimensions(prhs[0]);
  dL = mxGetNumberOfDimensions(prhs[0]);
  nL = mxGetNumberOfElements(prhs[0]);
  x  = (int)sL[0];
  y  = (int)sL[1];
  xy = x*y;

  sS[0] = 1; 
  sS[1] = 2; 
  SS = mxCreateNumericArray(2,sS,mxDOUBLE_CLASS,mxREAL);
  S = mxGetPr(SS);
  if (nrhs<3) {S[0]=1.0; S[1]=1.0; S[2]=1.0;} else {S=mxGetPr(prhs[2]);}
  
  s1 = fabs((float)S[0]),s2 = fabs((float)S[1]),s3 = fabs((float)S[2]);
  s12  = sqrt( s1*s1  + s2*s2); /* xy - voxel size */
  s13  = sqrt( s1*s1  + s3*s3); /* xz - voxel size */
  s23  = sqrt( s2*s2  + s3*s3); /* yz - voxel size */
  s123 = sqrt(s12*s12 + s3*s3); /* xyz - voxel size */
  /*printf("%1.2f,%1.2f,%1.2f - %1.2f,%1.2f,%1.2f - %1.2f",s1,s2,s3,s12,s23,s13,s123); */
  
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  NI[0] = 0;  NI[1] = -1;NI[2] = -x+1;NI[3] =  -x;NI[4] = -x-1;NI[5] =  -xy+1;NI[6] = -xy;NI[7] = -xy-1;NI[8] = -xy+x+1;NI[9] = -xy+x;NI[10] = -xy+x-1;NI[11] = -xy-x+1;NI[12] = -xy-x;NI[13] = -xy-x-1;  
  ND[0] = 0.0;ND[1] = s1;ND[2] =  s12;ND[3] =  s2;ND[4] =  s12;ND[5] =    s13;ND[6] = s3; ND[7] =   s13;ND[8] =    s123;ND[9] =   s23;ND[10] =    s123;ND[11] =    s123;ND[12] =   s23;ND[13] =   s123;
  sN = sizeof(NI)/4;    
  DN = (float *) mxMalloc(sizeof(float)*sN);
  DNm = FLT_MAX;
  i =0; n = 0; ni = 0; DNi = 0;

  /* data */
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[1] = mxCreateNumericArray(dL,sL,mxUINT32_CLASS,mxREAL);
  plhs[2] = mxCreateNumericArray(dL,sL,mxUINT8_CLASS,mxREAL);
  D = (float *)mxGetPr(plhs[0]);
  I = (unsigned int  *)mxGetPr(plhs[1]);
  L = (unsigned char *)mxGetPr(plhs[2]);
  
  V = (float*)mxGetPr(prhs[0]);
  if (nrhs>1) R=(bool *)mxGetPr(prhs[1]); 
  e255 = false; 

  /* intitialisation */
  for (i=0;i<nL;i++) 
  {
    if (V[i]>=0.5) D[i]=0.0; else D[i]=FLT_MAX; 
    if (V[i]>255.0)  
    {
      if (e255==false) 
      {
        printf("Warning: First parameter of vbdist > 255!\n"); 
        e255 = true;
      }
      V[i] = 255;
      
    }
    L[i]=(unsigned char) ceil(V[i]);
    I[i]=(unsigned int)i;
  }

  for (i=0;i<nL;i++) 
  {
    if ( (D[i]>0) && (nrhs==1 || (nrhs>1 && R[i]==true) ) )
    {
      ind2sub(i,&u,&v,&w,xy,x);
      
      /* read neighbor values */
      for (n=0;n<sN;n++)
      {
        ni = i + NI[n];
        ind2sub(ni,&nu,&nv,&nw,xy,x);
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) ) ni=i;
        DN[n] = D[ni] + ND[n];
      }

      /* find minimum distance within the neighborhood */
      pmin(DN,sN,&DNm,&DNi);

      /* update values */
      if (DNi>0) {
        L[i] = L[i+NI[DNi]];
        I[i] = (unsigned int) I[i+NI[DNi]];
        D[i] = DNm; 
        ind2sub((int)I[i],&nu,&nv,&nw,xy,x); 
        D[i] = sqrt(pow((float)(u-nu)*s1,2) + pow((float)(v-nv)*s2,2) + pow((float)(w-nw)*s3,2));
        /* hier muss die genauerer Berechnung mit pve rein! */
      }
    }
  }
  for (i=nL-1;i>0;i--)
  {
    if ( (D[i]>0) && (nrhs==1 || (nrhs>1 && R[i]==true) ) )
    {
      ind2sub(i,&u,&v,&w,xy,x);

      /* read neighbor values */
      for (n=0;n<sN;n++)
      {
        ni = i - NI[n];
        ind2sub(ni,&nu,&nv,&nw,xy,x);
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) ) ni=i;
        DN[n] = D[ni] + ND[n];
      }

      /* find minimum distance within the neighborhood */
      pmin(DN,sN,&DNm,&DNi);

      /* update values */
      if (DNi>0) {
        L[i] = L[i-NI[DNi]];
        I[i] = (unsigned int)  I[i-NI[DNi]];
        D[i] = DNm; 
        ind2sub((int)I[i],&nu,&nv,&nw,xy,x); 
        D[i] = sqrt(pow((float)(u-nu)*s1,2) + pow((float)(v-nv)*s2,2) + pow((float)(w-nw)*s3,2));
      }
    }
  }

  
  /* euclidean calcuation + PVE information */
  for (i=0;i<nL;i++) 
  {
  /*  if ( (D[i]>0) && (nrhs==1 || (R[i]>=1) ) )  D[i] = D[i] + (1 - V[I[i]]);
    if (D[i]<=0.5) D[i]=0; else D[i]=D[i]-0.5; */
    I[i]++;
  }

  mxFree(DN);
}


