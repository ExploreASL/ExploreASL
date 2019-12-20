/* Project-based volume (PBT) estimation 
 * _____________________________________________________________________________
 *
 * This function use the WM distance WMD to define a successor relation
 * ship in the cortex. 
 *
 *   [GMT,LV,RPM,LVc] = cat_vol_(SEG,WMD,CSFD)
 *
 *   SEG  = (single) segment image with low and high boundary bd
 *   WMD  = (single) CSF distance map
 *   CSFD = (single) CSF distance map
 *
 *   GMT  = (single) thickness image
 *   RPM  = (single) radial position map (equi-dist)
 * 
 * _____________________________________________________________________________
 * Robert Dahnke 201803
 * Center of Neuroimaging 
 * University Jena
 */

#include "mex.h"   
#include "matrix.h"
#include "math.h"
#include <stdlib.h>

struct opt_type {
	int   CSFD;													/* use CSFD */
	int   PVE;													/* 0, 1=fast, 2=exact */
	float LB, HB, LLB, HLB, LHB, HHB;  	/* boundary */
	int   sL[3];
	// ...
	} opt;

/* Function to estimate if neighbor *i are successor of *I given by the WM 
 * distance (WMDi and WMDI). Successors are neighbors of *I that have a  
 * specificly larger WM distance than the voxel itself, which is defined 
 * by the a1 and a2 boundary variables with 0 < a1 < 1 < a2 < 2. Small a1 
 * and large a2 increase the amount of voxels that are defined as successor
 * and create smoother results. 
 * Furthermore, it is expectd that the voxel *I and its successors *i should 
 * be in the GM ( projection range: round(SEGI)==2 and round(SEGi)==2 ). 
 * Moreover, we expect that the intensity of the successors *i are lighly
 * smaller that of *I. 
 */
bool issuccessor(float GMTi, float SEGi, float SEGI, 
                 float WMDi, float NDi, float WMDI,
                 float a1, float a2, float a3, float maximum) {
  return  ( ( GMTi < 1e15 ) && //( maximum < GMTi ) &&                /* thickness/WMD of neighbors should be larger */
            ( SEGi >  1.5 && SEGi <= 2.5 ) &&                         /* projection range - from  (1.2 - 2.75) */
            ( SEGI >  1.5 && SEGI <= 2.5 ) &&                         /* projection range - to    (1.5 - 2.75) */
            ( ( ( WMDi - NDi * a2 ) <= WMDI ) ) &&                    /* upper boundary - maximum distance */
            ( ( ( WMDi - NDi * a1 ) >  WMDI ) ) &&                    /* lower boundary - minimum distance - corrected values outside */
            ( ( ( SEGI - SEGi ) > (-a3 - (SEGi-1.5)/3 ) &&            /* SEGi should be not to much higher than SEGI */
                ( SEGI - SEGi ) < 0.5 ) || SEGI>2.0 )  ); 
}  
   

/* Find all successor voxels *i of a voxel *I to map the average thickness 
 * of the successors to *I.
 */
void pmax(const float GMT[], const float WMD[], const float SEG[], const float ND[], 
          const float WMDI,  const float CSFD,  const float SEGI, const int sA, float & maximum) {
  
  float T[27]; for (int i=0;i<27;i++) T[i]=-1; float n=0.0;
  maximum = WMDI; 
  
  float a1 =  0.5;   // lower boundary (lower values with include more voxels as successor) 0.6
  float a2 =  1.5;   // upper boundary (higher values with include more voxels as successor) 1.3
  float a3 =  0.15;  // minimum intensitiy difference between voxels (higher more successor)
  
  /* estiamte the maximum of sibblings */
  /* project volume values and count the siblings */
  for (int i=0;i<=sA;i++) {
    if ( issuccessor(GMT[i], SEG[i], SEGI, WMD[i], ND[i], WMDI, a1, a2, a3, maximum) ) {
      maximum = GMT[i];
    }
  }
  /* use the shortest distance between WM and CSF */
  maximum = fmin( maximum , WMDI + CSFD);
  
  /* the mean of the highest values of the siblings */
  float maximum2=maximum; float m2n=0.0; 
  for (int i=0;i<=sA;i++) {
    if ( issuccessor(GMT[i], SEG[i], SEGI, WMD[i], ND[i], WMDI, a1, a2, a3, maximum) ) {
      maximum2 = maximum2 + GMT[i]; m2n++; 
    }
  }
  if ( m2n > 0 )  maximum = (maximum2 - maximum) / m2n;

}




/* Estimate x,y,z position of index i in an array of size sx, sxy=sx*sy.
 * C index value 0,..,n-1 is used here rather than MATLAB index 1,..,n !
 */
void ind2sub(int i, int *x, int *y, int *z, int snL, int sxy, int sy) {
  /* not here ... 
   *  if (i<0) i=0; 
   *  if (i>=snL) i=snL-1;
  */
  
  *z = (int)floor( (double)i / (double)sxy ) ; 
   i = i % (sxy);
  *y = (int)floor( (double)i / (double)sy ) ;        
  *x = i % sy ;
}



/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs<3) mexErrMsgTxt("ERROR: not enought input elements\n");
  if (nrhs>4) mexErrMsgTxt("ERROR: to many input elements.\n");
  if (nlhs>4) mexErrMsgTxt("ERROR: to many output elements.\n");
  if (mxIsSingle(prhs[0])==0) mexErrMsgTxt("ERROR: first  input must be an 3d single matrix\n");
 
  
  /* main informations about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]); 
  mwSize sSEG[] = {sL[0],sL[1],sL[2]}; 
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     x  = sL[0];
  const int     y  = sL[1];
  const int     xy = x*y;
  const float   s2 = sqrt(2.0);
  const float   s3 = sqrt(3.0);
  const int     nr = nrhs;
  
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  const int   NI[]  = {  0, -1,-x+1, -x,-x-1,  -xy+1,-xy,-xy-1,  -xy+x+1,-xy+x,-xy+x-1,  -xy-x+1,-xy-x,-xy-x-1};  
  const float ND[]  = {0.0,1.0,  s2,1.0,  s2,     s2,1.0,   s2,       s3,   s2,     s3,       s3,   s2,     s3};
  const int   sN  = sizeof(NI)/4;  
  float       DN[sN],DI[sN],GMTN[sN],WMDN[sN],SEGN[sN],DNm,VOLN[sN],LVm;
  float*VOLc[sN]; 
  
  float 	    du, dv, dw, dnu, dnv, dnw, d, dcf, WMu, WMv, WMw, GMu, GMv, GMw, SEGl, SEGu, tmpfloat;
  int         mi,ni,u,v,w,nu,nv,nw, tmpint, WMC=0, CSFC=0;
    
  /* main volumes - actual without memory optimation ... */
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[1] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  
  /* input variables */
  float*SEG  = (float *)mxGetPr(prhs[0]);
  float*WMD  = (float *)mxGetPr(prhs[1]);
  float*CSFD = (float *)mxGetPr(prhs[2]);
  
  /* dynamic input did not work >> used hard code */
  /*if ( nrhs>1) {
		tmpint   = (int)mxGetScalar(mxGetField(prhs[1],1,"CSFD"));  printf("X=%d", tmpint);	if ( tmpint!=NULL && (tmpint>=0 && tmpint<=1) ) opt.CSFD = tmpint;   else opt.CSFD 	= 1;
		tmpint   = (int)mxGetScalar(mxGetField(prhs[1],1,"PVE"));   printf("X=%d", tmpint);	if ( tmpint!=NULL && (tmpint>=0 && tmpint<=2) ) opt.PVE  = tmpint; 	 else opt.PVE		= 2;
		tmpfloat = (float)mxGetScalar(mxGetField(prhs[1],1,"LB"));  printf("X=%d", tmpfloat);	if ( tmpfloat!=NULL ) 													opt.LB   = tmpfloat; else opt.LB  	= 1.5;
		tmpfloat = (float)mxGetScalar(mxGetField(prhs[1],1,"HB"));  printf("X=%d", tmpfloat);	if ( tmpfloat!=NULL ) 													opt.HB   = tmpfloat; else opt.HB 		= 2.5;
	} 
	else */{ opt.CSFD = 1;opt.PVE = 2;opt.LB = 1.5;opt.HB	= 2.5; }
	opt.LLB=floor(opt.LB), opt.HLB=ceil(opt.LB), opt.LHB=floor(opt.HB), opt.HHB=ceil(opt.HB);
  
  /* output variables */
  float        *GMT  = (float *)mxGetPr(plhs[0]);
  float        *RPM  = (float *)mxGetPr(plhs[1]);
  
  
  /* intitialisiation */
  for (int i=0;i<nL;i++) {
  	GMT[i] = WMD[i];
    RPM[i] = WMD[i];
		/* proof distance input */
    if ( SEG[i]>=opt.HB ) WMC++;
    if ( SEG[i]<=opt.LB ) CSFC++;
  }
	if (WMC==0)  mexErrMsgTxt("ERROR: no WM voxel\n");
	if (CSFC==0) opt.CSFD = 0;
  
  
/* Thickness mapping
 *  
 */
/* ============================================================================= */
  for (int i=0;i<nL;i++) {
    if (SEG[i]>opt.LLB && SEG[i]<opt.HHB) {
      ind2sub(i,&u,&v,&w,nL,xy,x);

      /* read neighbor values 
       * - why didn't you used pointers? 
       * - why not adding the neighbor loop to the subfunction?
       * > create an addition function 
       * > prepare convertation to C
       */
      for (int n=0;n<sN;n++) {
        ni = i + NI[n];
        ind2sub(ni,&nu,&nv,&nw,nL,xy,x);
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;
        GMTN[n] = GMT[ni]; WMDN[n] = RPM[ni]; SEGN[n] = SEG[ni]; 
      }

      /* find minimum distance within the neighborhood - forward */
      pmax(GMTN,WMDN,SEGN,ND,WMD[i],CSFD[i],SEG[i],sN,DNm);
      GMT[i] = DNm;
    }
  }

  for (int i=nL-1;i>=0;i--) {
    if (SEG[i]>opt.LLB && SEG[i]<opt.HHB) {
      ind2sub(i,&u,&v,&w,nL,xy,x);

      /* read neighbor values */
      for (int n=0;n<sN;n++) {
        ni = i - NI[n];
        ind2sub(ni,&nu,&nv,&nw,nL,xy,x);
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;
        GMTN[n] = GMT[ni]; WMDN[n] = RPM[ni]; SEGN[n] = SEG[ni]; 
      }

      /* find minimum distance within the neighborhood - backward */
      pmax(GMTN,WMDN,SEGN,ND,WMD[i],CSFD[i],SEG[i],sN,DNm);
      if ( GMT[i] < DNm && DNm>0 ) GMT[i] = DNm; 
    }
  }
  

  /* final GMT settings */
	for (int i=0;i<nL;i++) { 
    if (SEG[i]<opt.LB & SEG[i]>opt.LB)
      GMT[i] = 0;
    else
      GMT[i] = fmin( GMT[i], CSFD[i] + WMD[i] ); 
	}
 
  /* estimate RPM */
	for (int i=0;i<nL;i++) {
		if ( SEG[i]>=opt.HB ) 	
      RPM[i]=1.0; 
		else {
			if ( SEG[i]<=opt.LB || GMT[i]==0.0 ) 
        RPM[i]=0.0;
			else {
				RPM[i] = (GMT[i] - WMD[i]) / GMT[i];
				if (RPM[i]>1.0) RPM[i]=1.0;
				if (RPM[i]<0.0) RPM[i]=0.0;	
			}
		} 
	}
  
}



