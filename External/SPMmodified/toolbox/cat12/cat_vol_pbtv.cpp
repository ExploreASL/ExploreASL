/* Project-based volume (PBT) estiamtion 
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
 *   LV   = (single) lower volume map that descibe the "volumetric 
 *                   distance" rather than the distance
 *   LVc  = (single) count the number of successors of a voxel
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

/* some options */
struct opt_type {
	int   CSFD;													/* use CSFD */
	int   PVE;													/* 0, 1=fast, 2=exact */
	float LB, HB, LLB, HLB, LHB, HHB;  	/* boundary */
	int   sL[3];
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
bool issuccessor(float SEGi, float SEGI, float WMDi, float WMDI, float NDi, 
                 float a1, float a2, float a3) {
  return  ( ( SEGi >  1.5 && SEGi <= 2.5 ) &&                         /* projection range - from  (1.2 - 2.75) */
            ( SEGI >  1.5 && SEGI <= 2.5 ) &&                         /* projection range - to    (1.5 - 2.75) */
            ( ( ( WMDi - NDi * a2 ) <= WMDI ) ) &&                    /* upper boundary - maximum distance */
            ( ( ( WMDi - NDi * a1 ) >  WMDI ) ) &&                    /* lower boundary - minimum distance - corrected values outside */
            ( ( ( SEGI - SEGi ) < a3 ) || SEGI>1.8 ) );               /* SEGi should be not to much higher than SEGI */                  
            //( ( ( SEGI - SEGi ) > (-a3 - (SEGi-1.5)/3 ) &&            // SEGi should be not to much higher than SEGI 
}  
bool isprecessor(float SEGi, float SEGI, float WMDi, float WMDI, float NDi, 
                 float a1, float a2, float a3) {
  return  ( ( SEGi >  1.5 && SEGi <= 2.5 ) &&                         /* projection range - from  (1.2 - 2.75) */
            ( SEGI >  1.5 && SEGI <= 2.5 ) &&                         /* projection range - to    (1.5 - 2.75) */
            ( ( ( WMDI - NDi * a2 ) <= WMDi ) ) &&                    /* upper boundary - maximum distance */
            ( ( ( WMDI - NDi * a1 ) >  WMDi ) ) );                    /* lower boundary - minimum distance - corrected values outside */
  //          ( ( ( SEGi - SEGI ) < a3 ) || SEGi>1.8 ) );               /* SEGi should be not to much higher than SEGI */ 
  //( ( ( SEGi - SEGI ) > (-a3 - (SEGI-1.5)/3 ) &&     
}
bool isneighbor(float SEGi, float SEGI, float WMDi, float WMDI, float NDi, 
                float a1, float a2, float a3) {
  return  ( ( SEGi >  1.5 && SEGi <= 2.5 ) &&                         /* projection range - from  (1.2 - 2.75) */
            ( SEGI >  1.5 && SEGI <= 2.5 ) &&                         /* projection range - to    (1.5 - 2.75) */
            ( fabs( WMDi - WMDI ) <= a1  ) &&                         /* upper boundary - maximum distance */
            ( fabs( SEGi - SEGI ) <= a3  ) );                         /* SEGi should be not to much higher than SEGI */
}  

   

/* Find all successor voxels *i of a voxel *I to map the average thickness 
 * of the successors to *I.
 */
bool issuccessormax(float GMTi, float SEGi, float SEGI, 
                 float WMDi, float NDi, float WMDI,
                 float a1, float a2, float a3, float maximum) {
  return  ( ( GMTi < 1e15 ) && //( maximum < GMTi ) &&                /* thickness/WMD of neighbors should be larger */
            ( SEGi >  1.5 && SEGi <= 2.5 ) &&                         /* projection range - from  (1.2 - 2.75) */
            ( SEGI >  1.5 && SEGI <= 2.5 ) &&                         /* projection range - to    (1.5 - 2.75) */
            ( ( ( WMDi - NDi * a2 ) <= WMDI ) ) &&                    /* upper boundary - maximum distance */
            ( ( ( WMDi - NDi * a1 ) >  WMDI ) ) &&                    /* lower boundary - minimum distance - corrected values outside */
            ( ( ( SEGI - SEGi ) < a3 ) || SEGI>1.8 )  );              /* SEGi should be not to much higher than SEGI */
            //( ( ( SEGI - SEGi ) > (-a3 - (SEGi-1.5)/3 ) &&            /* SEGi should be not to much higher than SEGI */
}  
 

float pmax2(const float GMT[], const float WMD[], const float SEG[], const float ND[], 
          const float WMDI,  const float CSFD,  const float SEGI, const int sA) {
  
  float a1 =  0.5;   // lower boundary (lower values with include more voxels as successor) 0.6
  float a2 =  1.5;   // upper boundary (higher values with include more voxels as successor) 1.3
  float a3 =  0.1;  // minimum intensitiy difference between voxels (higher more successor)
  
  float m2n=1.0, maximum = WMDI; 
  for (int i=0;i<=sA;i++) {
    if ( (SEGI>SEG[i]+0.05) && issuccessor(SEG[i], SEGI, WMD[i], WMDI, ND[i], a1, a2, a3) ) {
      maximum = maximum + GMT[i]; m2n++; 
    }
  }
  return maximum / m2n;
}

/* Find all successor voxels *i of a voxel *I to map the average thickness 
 * of the successors to *I.
 */
float pmax(const float GMT[], const float WMD[], const float SEG[], const float ND[], 
          const float WMDI,  const float CSFD,  const float SEGI, const int sA) { //, float & maximum) {
  
  float T[27]; for (int i=0;i<27;i++) T[i]=-1; float n=0.0;
  float maximum = WMDI; 
  
  float a1 =  0.5;   // lower boundary (lower values with include more voxels as successor) 0.6
  float a2 =  1.5;   // upper boundary (higher values with include more voxels as successor) 1.3
  float a3 =  0.2;   // minimum intensitiy difference between voxels (higher more successor)
  
  /* estiamte the maximum of sibblings */
  /* project volume values and count the siblings */
  for (int i=0;i<=sA;i++) {
    if ( issuccessormax(GMT[i], SEG[i], SEGI, WMD[i], ND[i], WMDI, a1, a2, a3, maximum) ) {
      maximum = GMT[i];
    }
  }
  /* use the shortest distance between WM and CSF */
  maximum = fmin( maximum , WMDI + CSFD);
  
  /* the mean of the highest values of the siblings */
  float maximum2=maximum; float m2n=0.0; 
  for (int i=0;i<=sA;i++) {
    if ( issuccessormax(GMT[i], SEG[i], SEGI, WMD[i], ND[i], WMDI, a1, a2, a3, maximum) ) {
      maximum2 = maximum2 + GMT[i]; m2n++; 
    }
  }
  if ( m2n > 0 )  maximum = (maximum2 - maximum) / m2n;
  return maximum;
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
  if (nlhs>6) mexErrMsgTxt("ERROR: to many output elements.\n");
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
  plhs[2] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[3] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[4] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[5] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  
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
  float        *LV   = (float *)mxGetPr(plhs[1]); /* volume sum from CSF boundary */
  float        *UV   = (float *)mxGetPr(plhs[2]); /* volume sum from CSF boundary */
  float        *RPM  = (float *)mxGetPr(plhs[3]);
  float        *LVc  = (float *)mxGetPr(plhs[4]); /* volume sum counter */
  float        *UVc  = (float *)mxGetPr(plhs[5]); /* volume sum counter */
  
  
  /* intitialisiation */
  for (int i=0;i<nL;i++) {
  	GMT[i] = WMD[i];
    RPM[i] = WMD[i];
    /* for the LV map a partial volume initialization is used */
		LV[i]  = (float) fmin( 1 , fmin( CSFD[i] , WMD[i] )); // 1 - fmin( 1 , fabs( ( SEG[i]-2 ) * 2 ));
    UV[i]  = (float) fmin( 1 , fmin( CSFD[i] , WMD[i] )); // 1 - fmin( 1 , fabs( ( SEG[i]-2 ) * 2 ));
    LVc[i] = (float) (SEG[i]>=opt.LB && SEG[i]<=opt.HB); 
    UVc[i] = LVc[i]; 
		/* proof distance input */
    if ( SEG[i]>=opt.HB ) WMC++;
    if ( SEG[i]<=opt.LB ) CSFC++;
  }
	if (WMC==0)  mexErrMsgTxt("ERROR: no WM voxel\n");
	if (CSFC==0) opt.CSFD = 0;
  
  
  
  /* estimate successor and precessor weighting */
  float a1 =  0.1;   // lower boundary (lower values with include more voxels as successor) 0.6
  float a2 =  2.5;   // upper boundary (higher values with include more voxels as successor) 1.3
  float a3 =  0.5;   // minimum intensitiy difference between voxels (higher more successor)
  for (int i=0;i<nL;i++) {
    if (SEG[i]>opt.LLB && SEG[i]<opt.HHB) {
      ind2sub(i,&u,&v,&w,nL,xy,x);

      for (int di=0;di<2;di++) { 
        for (int n=0;n<sN;n++) {
          if ( di==0 ) ni = i + NI[n]; else ni = i - NI[n]; 
          ind2sub(ni,&nu,&nv,&nw,nL,xy,x);
          if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;

          LVc[i] = LVc[i] + (float)(issuccessor(SEG[ni], SEG[i], WMD[ni], WMD[i], ND[n], a1, a2, a3));
          UVc[i] = UVc[i] + (float)(isprecessor(SEG[ni], SEG[i], WMD[ni], WMD[i], ND[n], a1, a2, a3));
        }
      }
    }
  }
          

   
  /* volume mapping - CSF to WM */
  int ip; 
  for (int di=0;di<2;di++) {
    if (di==0) ip = 1; else ip = -1;
    for (int i=di*(nL-1);i>=0 && i<nL;i=i+ip) {
      if (SEG[i]>opt.LLB && SEG[i]<opt.HHB) {
        ind2sub(i,&u,&v,&w,nL,xy,x);
        for (int n=0;n<sN;n++) {
          if ( di==0 ) ni = i + NI[n]; else ni = i - NI[n]; 
          ind2sub(ni,&nu,&nv,&nw,nL,xy,x);
          if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;
         
          if ( issuccessor(SEG[ni], SEG[i], WMD[ni], WMD[i], ND[n], a1, a2, a3) ) 
            LV[i] = LV[i] + LV[ni] / fmax(1,UVc[ni]); 
          if ( isprecessor(SEG[ni], SEG[i], WMD[ni], WMD[i], ND[n], a1, a2, a3) )
            UV[i] = UV[i] + UV[ni] / fmax(1,LVc[ni]);
        }
      }
    }
  }
  /* smoothing */
  float nlc, nuc, rpd;  
  for (int dic=0;dic<4;dic++) {
    for (int di=0;di<2;di++) {
      if (di==0) ip = 1; else ip = -1;
      for (int i=di*(nL-1);i>=0 && i<nL;i=i+ip) {
        if (SEG[i]>opt.LLB && SEG[i]<opt.HHB) {
          ind2sub(i,&u,&v,&w,nL,xy,x);
          nlc=1.0; nuc=1.0;  
          for (int n=0;n<sN;n++) {
            if ( di==0 ) ni = i + NI[n]; else ni = i - NI[n]; 
            ind2sub(ni,&nu,&nv,&nw,nL,xy,x);
            if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;

            if ( isneighbor(SEG[ni], SEG[i], WMD[ni], WMD[i], ND[n], 0.4, a2, 0.2) && 
                 (LV[i]/(LV[i] + UV[i]))>0.1 && (UV[i]/(LV[i] + UV[i]))>0.1  ) {
            
              rpd  = fmax(0,1 - fabs( LV[i]/(LV[i] + UV[i]) - LV[ni]/(LV[ni] + UV[ni]) * 4))*0.5; 
              rpd += fmax(0,1 - fabs( SEG[i] - SEG[ni] ) * 10)*0.5; 
              
              if (rpd>0) {
                LV[i] = LV[i] + LV[ni]*rpd; 
                nlc += rpd;
              }
              rpd  = fmax(0,1 - fabs( UV[i]/(LV[i] + UV[i]) - UV[ni]/(LV[ni] + UV[ni]) * 4))*0.5;
              rpd += fmax(0,1 - fabs( SEG[i] - SEG[ni] ) * 10)*0.5; 
              if ( rpd>0 ) {
                UV[i] = UV[i] + UV[ni]*rpd; 
                nuc += rpd;
              }
            }
          }
          LV[i] = LV[i] / nlc;
          UV[i] = UV[i] / nuc;
        }
      }
    }
  }
  
  

  /* WMD correction */
  for (int i=0;i<nL;i++) {
    if ( ( LVc[i] > 2.0 ) && ( UVc[i] < 20.0 ) && ( LVc[i] > UVc[i] ) && 
         ( LV[i]  > 6.0 ) && ( UV[i]  <  2.0 ) && ( LV[i]  > UV[i]  ) && 
         ( SEG[i] > 1.9 ) && ( SEG[i] <  2.5 ) )   
      WMD[i] = 0.5; 
  }
  for (int di=0;di<2;di++) {
    if (di==0) ip = 1; else ip = -1;
    for (int i=di*(nL-1);i>=0 && i<nL;i=i+ip) {
      if ( WMD[i]>0 && SEG[i]<2.5) {
        ind2sub(i,&u,&v,&w,nL,xy,x);

        /* read neighbor values 
         * - why didn't you used pointers? 
         * - why not adding the neighbor loop to the subfunction?
         * > create an addition function 
         * > prepare convertation to C
         */
        for (int n=0;n<sN;n++) {
          if ( di==0 ) ni = i + NI[n]; else ni = i - NI[n]; 
          ind2sub(ni,&nu,&nv,&nw,nL,xy,x);
          if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;
          if ( ( WMD[ni]>=0.4 ) && ( WMD[i]>=1 ) && ( WMD[i] > WMD[ni] + ND[n]) && 
               ( SEG[ni]<2.5 ) && ( SEG[ni]>1.5 ) && ( SEG[i]<=SEG[ni]+0.2 ) )
            WMD[i] = WMD[ni] + ND[n]; 
        }
      }
    }
  }
  
  
  
  // printf("Thickness mapping\n");
  /* Thickness mapping */
  /* ============================================================================= */
if ( 0 ) {
  /* new shorter version that did not run correctly*/
  for (int di=0;di<2;di++) {
    if (di==0) ip = 1; else ip = -1;
    for (int i=di*(nL-1);i>=0 && i<nL;i=i+ip) {
      if (SEG[i]>opt.LLB && SEG[i]<opt.HHB) {
        ind2sub(i,&u,&v,&w,nL,xy,x);

        /* read neighbor values 
         * - why didn't you used pointers? 
         * - why not adding the neighbor loop to the subfunction?
         * > create an addition function 
         * > prepare convertation to C
         */
        for (int n=0;n<sN;n++) {
          if ( di==0 ) ni = i + NI[n]; else ni = i - NI[n]; 
          ind2sub(ni,&nu,&nv,&nw,nL,xy,x);
          if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;
          GMTN[n] = GMT[ni]; WMDN[n] = RPM[ni]; SEGN[n] = SEG[ni]; 
        }

        // find minimum distance within the neighborhood - forward 
        if ( di == 0)
          GMT[i] = pmax2(GMTN,WMDN,SEGN,ND,WMD[i],CSFD[i],SEG[i],sN);
        else {
          DNm = pmax2(GMTN,WMDN,SEGN,ND,WMD[i],CSFD[i],SEG[i],sN);
          if ( GMT[i] < DNm && DNm>0 ) GMT[i] = DNm; 
         }
      }
    }
  }
}
else {
  /* old long code version that works */
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
      
      GMT[i] = pmax(GMTN,WMDN,SEGN,ND,WMD[i],CSFD[i],SEG[i],sN);
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
      DNm = pmax(GMTN,WMDN,SEGN,ND,WMD[i],CSFD[i],SEG[i],sN);
      if ( GMT[i] < DNm && DNm>0 ) GMT[i] = DNm; 
    }
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
		if ( SEG[i]>opt.HB ) 	
      RPM[i]=1.0; 
		else {
			if ( SEG[i]<opt.LB || GMT[i]==0.0 ) 
        RPM[i]=0.0;
			else {
				RPM[i] = (GMT[i] - WMD[i]) / GMT[i];
				if (RPM[i]>1.0) RPM[i]=1.0;
				if (RPM[i]<0.0) RPM[i]=0.0;	
			}
		} 
	}
  
}


