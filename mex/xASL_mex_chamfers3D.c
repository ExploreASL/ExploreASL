/*
 * xASL_mex_chamfers3D.c
 * Borgefors Chamfers computation of Euclidean distance in 3D using a 5x5x5 window
 * Binary image is on the input
 * Returned is the distance map, and three components (X,Y,Z) of the nearest voxel
 *
 * [dist,x,y,z] = xASL_mex_chamfers3D(im)
 * im - 3D image, double, binary mask
 * dist - 3D image of distances from the nearest voxel on the mask
 * x,y,z - 3D images - coordinates of the nearest voxel on the mask
 * 
 * Implemented by Jan Petr
 * According to: Stina Svensson and Gunilla Borgefors, Computer Vision and Image Understanding 88, 24–53 (2002)
 * doi:10.1006/cviu.2002.0976
 */

#include <string.h>
#include "mex.h"

#ifndef CDA
#define CDA 1
#endif

#ifndef CDB
#define CDB 1.4142
#endif

#ifndef CDC
#define CDC 1.7321
#endif

#ifndef CDD
#define CDD 2.2361
#endif

#ifndef CDE
#define CDE 2.4495
#endif

#ifndef CDF
#define CDF 3
#endif


/*#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif*/

/* Utility function that returns index into */
/* 1D array with range checking.            */
 
/* This function checks if the maximum distance of the current coordinate (X,Y,Z)
 * Cannot be minimized from the neighbor(X+I,Y+J,Z+K)
 * If yes, then the distance and the 3D vector to the nearest structure are updated*/

void update_dist(double *poimd,
            double *poimx,
            double *poimy,        
            double *poimz,
            mwSignedIndex i,
            mwSignedIndex j,
            mwSignedIndex k,
            mwIndex inei,
            double D)
{
    *(poimd) = *(poimd+inei) + D;
    *(poimx) = *(poimx+inei) + i;
    *(poimy) = *(poimy+inei) + j;
    *(poimz) = *(poimz+inei) + k;
}

void check_dist(double *poimd,
                double *poimx,
                double *poimy,        
                double *poimz,
                mwSignedIndex i,
                mwSignedIndex j,
                mwSignedIndex k,
                mwSize dim[3],
                double D)
{
    mwIndex inei;
    
    // Calculate the neighbor index and proceed only if within the image
    inei = i + j*dim[0] + k*dim[0]*dim[1];
    
    if ( *(poimd) > (*(poimd+inei) + D) )
        update_dist(poimd,poimx,poimy,poimz,i,j,k,inei,D);
}

// Function doing the job. 

void do_it(double  *iima,
           mwSize  dim[3],
           double  *oimd,
           double  *oimx,
           double  *oimy,
           double  *oimz)
{
    int x,y,z;
    int inei;
    int n,nx,nxy,n2x,n2xy;
    double *piima;
    double *poimd,*poimx,*poimy,*poimz;
    int z2,z1,zm1,zm2,x2,x1,xm1,xm2,y2,y1,ym1,ym2;
    int m1mn2xmn2xy,p1mn2xmn2xy,m2mnxmn2xy,
            m1mnxmn2xy,
            mnxmn2xy,
            p1mnxmn2xy,
            p2mnxmn2xy,
            m1mn2xy,
            p1mn2xy,
            m2pnxmn2xy,
            m1pnxmn2xy,
            pnxmn2xy,
            p1pnxmn2xy,
            p2pnxmn2xy,
            m1pn2xmn2xy,
            p1pn2xmn2xy,
            m2mn2xmnxy,
            m1mn2xmnxy,
            mn2xmnxy,
            p1mn2xmnxy,
            p2mn2xmnxy,
            m2mnxmnxy,
            m1mnxmnxy,
            mnxmnxy,
            p1mnxmnxy,
            p2mnxmnxy,
            m2mnxy,
            m1mnxy,
            p1mnxy,
            p2mnxy,
            m2pnxmnxy,
            m1pnxmnxy,
            pnxmnxy,
            p1pnxmnxy,
            p2pnxmnxy,
            m2pn2xmnxy,
            m1pn2xmnxy,
            pn2xmnxy,
            p1pn2xmnxy,
            p2pn2xmnxy,
            m1mn2x,
            p1mn2x,
            m2mnx,
            m1mnx,
            p1mnx,
            p2mnx;
    
    n = dim[0]*dim[1]*dim[2];
    nx = dim[0];
    nxy = dim[0]*dim[1];
    n2x = 2*dim[0];
    n2xy = 2*dim[0]*dim[1];
    m1mn2xmn2xy = -1-n2x-n2xy;
    p1mn2xmn2xy = 1-n2x-n2xy;
    m2mnxmn2xy = -2-nx-n2xy;
    m1mnxmn2xy = -1-nx-n2xy;
    mnxmn2xy = -nx-n2xy;
    p1mnxmn2xy = 1-nx-n2xy;
    p2mnxmn2xy = 2-nx-n2xy;
    m1mn2xy = -1-n2xy;
    p1mn2xy = 1-n2xy;
    m2pnxmn2xy = -2+nx-n2xy;
    m1pnxmn2xy = -1+nx-n2xy;
    pnxmn2xy = nx-n2xy;
    p1pnxmn2xy = 1+nx-n2xy;
    p2pnxmn2xy = 2+nx-n2xy;
    m1pn2xmn2xy = -1+n2x-n2xy;
    p1pn2xmn2xy = 1+n2x-n2xy;
    m2mn2xmnxy = -2-n2x-nxy;
    m1mn2xmnxy = -1-n2x-nxy;
    mn2xmnxy = -n2x-nxy;
    p1mn2xmnxy = 1-n2x-nxy;
    p2mn2xmnxy = 2-n2x-nxy;
    m2mnxmnxy = -2-nx-nxy;
    m1mnxmnxy = -1-nx-nxy;
    mnxmnxy = -nx-nxy;
    p1mnxmnxy = 1-nx-nxy;
    p2mnxmnxy = 2-nx-nxy;
    m2mnxy = -2-nxy;
    m1mnxy = -1-nxy;
    p1mnxy = 1-nxy;
    p2mnxy = 2-nxy;
    m2pnxmnxy = -2+nx-nxy;
    m1pnxmnxy = -1+nx-nxy;
    pnxmnxy = nx-nxy;
    p1pnxmnxy = 1+nx-nxy;
    p2pnxmnxy = 2+nx-nxy;
    m2pn2xmnxy = -2+n2x-nxy;
    m1pn2xmnxy = -1+n2x-nxy;
    pn2xmnxy = n2x-nxy;
    p1pn2xmnxy = 1+n2x-nxy;
    p2pn2xmnxy = 2+n2x-nxy;
    m1mn2x = -1-n2x;
    p1mn2x = 1-n2x;
    m2mnx = -2-nx;
    m1mnx = -1-nx;
    p1mnx = 1-nx;
    p2mnx = 2-nx;
    
	
    //mexPrintf("Size: %d,%d,%d\n",(int) dim[0],(int) dim[1],(int) dim[2]);
    
    // Initialize the distance and shift components to inf/zero and zero 
    piima = iima;
    poimd = oimd;
    poimx = oimx;
    poimy = oimy;
    poimz = oimz;
    for (x=0;x<n;x++)
    {
        *(poimx++) = 0;
        *(poimy++) = 0;
        *(poimz++) = 0;
        if (*piima++ > 0)
            *poimd = 0;
        else
            *poimd = n;
        poimd++;
    }
    
    // Forward loop 
    poimd = oimd;
    poimx = oimx;
    poimy = oimy;
    poimz = oimz;
    for (z=0;z<dim[2];z++)
    {
        for (y=0;y<dim[1];y++)
        {
            for (x=0;x<dim[0];x++)
            {
				
                // Only execute if the distance is greater than 1 as 0 or 1 cannot be improved 
                
                if ((*poimd) > 1)
                {
					
                    zm2 = (z>=2);
                    ym2 = (y>=2);
                    xm2 = (x>=2);
                    zm1 = (z>=1);
                    ym1 = (y>=1);
                    xm1 = (x>=1);
                    y1 = ((y+1)<dim[1]);
                    y2 = ((y+2)<dim[1]);
                    x1 = ((x+1)<dim[0]);
                    x2 = ((x+2)<dim[0]);
                    
                    //check_dist(poimd,poimx,poimy,poimz,i,j,k,dim,D);

                    if (zm2)
                    {
                        if (ym2)
                        {
                            if (xm1) {inei = m1mn2xmn2xy;if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz,-1,-2,-2,inei,CDF);}
                            if (x1)  {inei = p1mn2xmn2xy;if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz, 1,-2,-2,inei,CDF);}
                        }
                        if (ym1)
                        {
                            if (xm2) {inei = m2mnxmn2xy; if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz,-2,-1,-2,inei,CDF);}
                            if (xm1) {inei = m1mnxmn2xy; if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz,-1,-1,-2,inei,CDE);}
                                      inei =   mnxmn2xy; if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz, 0,-1,-2,inei,CDD);
                            if (x1)  {inei = p1mnxmn2xy; if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz, 1,-1,-2,inei,CDE);}
                            if (x2)  {inei = p2mnxmn2xy; if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz, 2,-1,-2,inei,CDF);}
                         
                        }
                            if (xm1) {inei =    m1mn2xy; if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz,-1,0,-2,inei,CDD);}
                            if (x1)  {inei =    p1mn2xy; if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz, 1,0,-2,inei,CDD);}
                        if (y1)
                        {
                            if (xm2) {inei = m2pnxmn2xy; if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz,-2,1,-2,inei,CDF);}
                            if (xm1) {inei = m1pnxmn2xy; if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz,-1,1,-2,inei,CDE);}
                                      inei =   pnxmn2xy; if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz, 0,1,-2,inei,CDD);
                            if (x1)  {inei = p1pnxmn2xy; if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz, 1,1,-2,inei,CDE);}
                            if (x2)  {inei = p2pnxmn2xy; if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz, 2,1,-2,inei,CDF);}
                        }
                        if (y2)
                        {
                            if (xm1) {inei = p1pn2xmn2xy;if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz,-1,2,-2,inei,CDF);}
                            if (x1)  {inei = p1pn2xmn2xy;if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz, 1,2,-2,inei,CDF);}
                        }
                    }
					
                    if (zm1)
                    {
                        if (ym2)
                        {
                            if (xm2) {inei =  m2mn2xmnxy; if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz,-2,-2,-1,inei,CDF);}
                            if (xm1) {inei =  m1mn2xmnxy; if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz,-1,-2,-1,inei,CDE);}
                                      inei =    mn2xmnxy; if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz, 0,-2,-1,inei,CDD);
                            if  (x1) {inei =  p1mn2xmnxy; if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz, 1,-2,-1,inei,CDE);}
                            if  (x2) {inei =  p2mn2xmnxy; if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz, 2,-2,-1,inei,CDF);}
                        }
                        if (ym1)
                        {
                            if (xm2) {inei =   m2mnxmnxy; if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz,-2,-1,-1,inei,CDE);}
                            if (xm1) {inei =   m1mnxmnxy; if (*(poimd)>(*(poimd+inei)+CDC)) update_dist(poimd,poimx,poimy,poimz,-1,-1,-1,inei,CDC);}
                                      inei =     mnxmnxy; if (*(poimd)>(*(poimd+inei)+CDB)) update_dist(poimd,poimx,poimy,poimz, 0,-1,-1,inei,CDB);
                            if (x1)  {inei =   p1mnxmnxy; if (*(poimd)>(*(poimd+inei)+CDC)) update_dist(poimd,poimx,poimy,poimz, 1,-1,-1,inei,CDC);}
                            if (x2)  {inei =   p2mnxmnxy; if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz, 2,-1,-1,inei,CDE);}
                        }
                        
                            if (xm2) {inei =      m2mnxy; if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz,-2, 0,-1,inei,CDD);}
                            if (xm1) {inei =      m1mnxy; if (*(poimd)>(*(poimd+inei)+CDB)) update_dist(poimd,poimx,poimy,poimz,-1, 0,-1,inei,CDB);}
                                      inei =        -nxy; if (*(poimd)>(*(poimd+inei)+CDA)) update_dist(poimd,poimx,poimy,poimz, 0, 0,-1,inei,CDA);
                            if (x1)  {inei =      p1mnxy; if (*(poimd)>(*(poimd+inei)+CDB)) update_dist(poimd,poimx,poimy,poimz, 1, 0,-1,inei,CDB);}
                            if (x2)  {inei =      p2mnxy; if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz, 2, 0,-1,inei,CDD);}
                        
                        if (y1)
                        {
                            if (xm2) {inei =   m2pnxmnxy; if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz,-2, 1,-1,inei,CDE);}
                            if (xm1) {inei =   m1pnxmnxy; if (*(poimd)>(*(poimd+inei)+CDC)) update_dist(poimd,poimx,poimy,poimz,-1, 1,-1,inei,CDC);}
                                      inei =     pnxmnxy; if (*(poimd)>(*(poimd+inei)+CDB)) update_dist(poimd,poimx,poimy,poimz, 0, 1,-1,inei,CDB);
                            if (x1)  {inei =   p1pnxmnxy; if (*(poimd)>(*(poimd+inei)+CDC)) update_dist(poimd,poimx,poimy,poimz, 1, 1,-1,inei,CDC);}
                            if (x2)  {inei =   p2pnxmnxy; if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz, 2, 1,-1,inei,CDE);}
                        }
                        if (y2)
                        {
                            if (xm2) {inei =  m2pn2xmnxy; if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz,-2, 2,-1,inei,CDF);}
                            if (xm1) {inei =  m1pn2xmnxy; if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz,-1, 2,-1,inei,CDE);}
                                      inei =    pn2xmnxy; if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz, 0, 2,-1,inei,CDD);
                            if (x1)  {inei =  p1pn2xmnxy; if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz, 1, 2,-1,inei,CDE);}
                            if (x2)  {inei =  p2pn2xmnxy; if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz, 2, 2,-1,inei,CDF);}
                        }
                    }
                    
                    if (ym2)
                    {
                            if (xm1) {inei =      m1mn2x; if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz,-1,-2, 0,inei,CDD);}
                            if (x1)  {inei =      p1mn2x; if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz, 1,-2, 0,inei,CDD);}
                    }
                    if (ym1)
                    {
                            if (xm2) {inei =       m2mnx; if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz,-2,-1, 0,inei,CDD);}
                            if (xm1) {inei =       m1mnx; if (*(poimd)>(*(poimd+inei)+CDB)) update_dist(poimd,poimx,poimy,poimz,-1,-1, 0,inei,CDB);}
                                      inei =         -nx; if (*(poimd)>(*(poimd+inei)+CDA)) update_dist(poimd,poimx,poimy,poimz, 0,-1, 0,inei,CDA);
                            if (x1)  {inei =       p1mnx; if (*(poimd)>(*(poimd+inei)+CDB)) update_dist(poimd,poimx,poimy,poimz, 1,-1, 0,inei,CDB);}
                            if (x2)  {inei =       p2mnx; if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz, 2,-1, 0,inei,CDD);}
                    }
                            if (xm1) {inei =          -1; if (*(poimd)>(*(poimd+inei)+CDA)) update_dist(poimd,poimx,poimy,poimz,-1, 0, 0,inei,CDA);}

                }
                poimd++;
                poimx++;
                poimy++;
                poimz++;
				
            }
        }
    }
    
    // Backward loop 
    poimd = oimd+n-1;
    poimx = oimx+n-1;
    poimy = oimy+n-1;
    poimz = oimz+n-1;
	
    for (z=(dim[2]-1);z>=0;z--)
    {
        for (y=dim[1]-1;y>=0;y--)
        {
            for (x=dim[0]-1;x>=0;x--)
            {
                // Only execute if the distance is greater than 1 as 0 or 1 cannot be improved 
                
                if ((*poimd) > 1)
                {
                    ym2 = (y>=2);
                    xm2 = (x>=2);
                    ym1 = (y>=1);
                    xm1 = (x>=1);
                    z1 = ((z+1)<dim[2]);
                    z2 = ((z+2)<dim[2]);
                    y1 = ((y+1)<dim[1]);
                    y2 = ((y+2)<dim[1]);
                    x1 = ((x+1)<dim[0]);
                    x2 = ((x+2)<dim[0]);
                    //check_dist(poimd,poimx,poimy,poimz,i,j,k,dim,D);

					if (z2)
                    {
                        if (ym2)
                        {
                            if (xm1) {inei = -1-n2x+n2xy;if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz,-1,-2, 2,inei,CDF);}
                            if (x1)  {inei =  1-n2x+n2xy;if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz, 1,-2, 2,inei,CDF);}
                        }
                        
                        if (ym1)
                        {
                            if (xm2) {inei =  -2-nx+n2xy;if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz,-2,-1, 2,inei,CDF);}
                            if (xm1) {inei =  -1-nx+n2xy;if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz,-1,-1, 2,inei,CDE);}
                                      inei =    -nx+n2xy;if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz, 0,-1, 2,inei,CDD);
                            if (x1)  {inei =   1-nx+n2xy;if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz, 1,-1, 2,inei,CDE);}
                            if (x2)  {inei =   2-nx+n2xy;if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz, 2,-1, 2,inei,CDF);}
                         
                        }
                            if (xm1) {inei =     -1+n2xy;if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz,-1, 0, 2,inei,CDD);}
                            if (x1)  {inei =      1+n2xy;if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz, 1, 0, 2,inei,CDD);}
                        if (y1)
                        {
                            if (xm2) {inei =  -2+nx+n2xy;if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz,-2, 1, 2,inei,CDF);}
                            if (xm1) {inei =  -1+nx+n2xy;if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz,-1, 1, 2,inei,CDE);}
                                      inei =     nx+n2xy;if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz, 0, 1, 2,inei,CDD);
                            if (x1)  {inei =   1+nx+n2xy;if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz, 1, 1, 2,inei,CDE);}
                            if (x2)  {inei =   2+nx+n2xy;if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz, 2, 1, 2,inei,CDF);}
                        }
                        if (y2)
                        {
                            if (xm1) {inei = -1+n2x+n2xy;if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz,-1, 2, 2,inei,CDF);}
                            if (x1)  {inei =  1+n2x+n2xy;if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz, 1, 2, 2,inei,CDF);}
                        }
                    }
                    
                    if (z1)
                    {
                        if (ym2)
                        {
                            if (xm2) {inei =  -2-n2x+nxy;if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz,-2,-2, 1,inei,CDF);}
                            if (xm1) {inei =  -1-n2x+nxy;if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz,-1,-2, 1,inei,CDE);}
                                      inei =    -n2x+nxy;if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz, 0,-2, 1,inei,CDD);
                            if (x1)  {inei =   1-n2x+nxy;if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz, 1,-2, 1,inei,CDE);}
                            if (x2)  {inei =   2-n2x+nxy;if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz, 2,-2, 1,inei,CDF);}
                        }
                        if (ym1)
                        {
                            if (xm2) {inei =   -2-nx+nxy;if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz,-2,-1, 1,inei,CDE);}
                            if (xm1) {inei =   -1-nx+nxy;if (*(poimd)>(*(poimd+inei)+CDC)) update_dist(poimd,poimx,poimy,poimz,-1,-1, 1,inei,CDC);}
                                      inei =     -nx+nxy;if (*(poimd)>(*(poimd+inei)+CDB)) update_dist(poimd,poimx,poimy,poimz, 0,-1, 1,inei,CDB);
                            if (x1)  {inei =    1-nx+nxy;if (*(poimd)>(*(poimd+inei)+CDC)) update_dist(poimd,poimx,poimy,poimz, 1,-1, 1,inei,CDC);}
                            if (x2)  {inei =    2-nx+nxy;if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz, 2,-1, 1,inei,CDE);}
                        }
                        
                            if (xm2) {inei =      -2+nxy;if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz,-2, 0, 1,inei,CDD);}
                            if (xm1) {inei =      -1+nxy;if (*(poimd)>(*(poimd+inei)+CDB)) update_dist(poimd,poimx,poimy,poimz,-1, 0, 1,inei,CDB);}
                                      inei =         nxy;if (*(poimd)>(*(poimd+inei)+CDA)) update_dist(poimd,poimx,poimy,poimz, 0, 0, 1,inei,CDA);
                            if (x1)  {inei =       1+nxy;if (*(poimd)>(*(poimd+inei)+CDB)) update_dist(poimd,poimx,poimy,poimz, 1, 0, 1,inei,CDB);}
                            if (x2)  {inei =       2+nxy;if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz, 2, 0, 1,inei,CDD);}
                        
                        if (y1)
                        {
                            if (xm2) {inei =   -2+nx+nxy;if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz,-2,1,1,inei,CDE);}
                            if (xm1) {inei =   -1+nx+nxy;if (*(poimd)>(*(poimd+inei)+CDC)) update_dist(poimd,poimx,poimy,poimz,-1,1,1,inei,CDC);}
                                      inei =      nx+nxy;if (*(poimd)>(*(poimd+inei)+CDB)) update_dist(poimd,poimx,poimy,poimz,0,1,1,inei,CDB);
                            if (x1)  {inei =    1+nx+nxy;if (*(poimd)>(*(poimd+inei)+CDC)) update_dist(poimd,poimx,poimy,poimz,1,1,1,inei,CDC);}
                            if (x2)  {inei =    2+nx+nxy;if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz,2,1,1,inei,CDE);}
                        }
                        if (y2)
                        {
                            if (xm2) {inei =  -2+n2x+nxy;if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz,-2, 2, 1,inei,CDF);}
                            if (xm1) {inei =  -1+n2x+nxy;if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz,-1, 2, 1,inei,CDE);}
                                      inei =     n2x+nxy;if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz, 0, 2, 1,inei,CDD);
                            if (x1)  {inei =   1+n2x+nxy;if (*(poimd)>(*(poimd+inei)+CDE)) update_dist(poimd,poimx,poimy,poimz, 1, 2, 1,inei,CDE);}
                            if (x2)  {inei =   2+n2x+nxy;if (*(poimd)>(*(poimd+inei)+CDF)) update_dist(poimd,poimx,poimy,poimz, 2, 2, 1,inei,CDF);}
                        }
                    }
                   
                    if (y2)
                    {
                            if (xm1) {inei =      -1+n2x;if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz,-1, 2, 0,inei,CDD);}
                            if (x1)  {inei =       1+n2x;if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz, 1, 2, 0,inei,CDD);}
                    }
                    if (y1)
                    {
                            if (xm2) {inei =       -2+nx;if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz,-2, 1, 0,inei,CDD);}
                            if (xm1) {inei =       -1+nx;if (*(poimd)>(*(poimd+inei)+CDB)) update_dist(poimd,poimx,poimy,poimz,-1, 1, 0,inei,CDB);}
                                      inei =          nx;if (*(poimd)>(*(poimd+inei)+CDA)) update_dist(poimd,poimx,poimy,poimz, 0, 1, 0,inei,CDA);
                            if (x1)  {inei =        1+nx;if (*(poimd)>(*(poimd+inei)+CDB)) update_dist(poimd,poimx,poimy,poimz, 1, 1, 0,inei,CDB);}
                            if (x2)  {inei =        2+nx;if (*(poimd)>(*(poimd+inei)+CDD)) update_dist(poimd,poimx,poimy,poimz, 2, 1, 0,inei,CDD);}
                    }
                            if (x1)  {inei =           1;if (*(poimd)>(*(poimd+inei)+CDA)) update_dist(poimd,poimx,poimy,poimz, 1, 0, 0,inei,CDA);}

                }	 
                poimd--;
                poimx--;
                poimy--;
                poimz--;
		
            }
        }
    } 
}
    
   


/* Gateway function */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize        ndim=0;
    int           i;
    const mwSize  *cdim = NULL;
    mwSize        dim[3];
    double        *iima = NULL;
    double        *oimd = NULL, *oimx = NULL, *oimy = NULL, *oimz = NULL;

    if (nrhs < 1) mexErrMsgTxt("xASL_mex_chamfers3D: Not enough input arguments.");
    if (nrhs > 1) mexErrMsgTxt("xASL_mex_chamfers3D: Too many input arguments.");
    if (nlhs < 4) mexErrMsgTxt("xASL_mex_chamfers3D: Not enough output arguments");
    if (nlhs > 4) mexErrMsgTxt("xASL_mex_chamfers3D: Too many output arguments.");

    /* Get image */

    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
    {
        mexErrMsgTxt("xASL_mex_chamfers3D: ima must be numeric, real, full and double");
    }
    ndim = mxGetNumberOfDimensions(prhs[0]);
    if ((ndim < 2) || (ndim > 3))
    {
        mexErrMsgTxt("xASL_mex_chamfers3D: ima must be 2- or 3-dimensional");
    }
    cdim = mxGetDimensions(prhs[0]);
    iima = mxGetPr(prhs[0]);

    
    /* Fix dimensions to allow for 2D and 3D data */

    dim[0] = cdim[0]; dim[1] = cdim[1];
    if (ndim==2) {dim[2]=1; ndim=3;} else {dim[2]=cdim[2];} 
  
    /* Allocate and initialise output images */

    plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
                                   mxGetDimensions(prhs[0]),mxDOUBLE_CLASS,mxREAL);
    oimd = mxGetPr(plhs[0]);
    
    plhs[2] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
                                   mxGetDimensions(prhs[0]),mxDOUBLE_CLASS,mxREAL);
    oimx = mxGetPr(plhs[2]);
    
    plhs[1] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
                                   mxGetDimensions(prhs[0]),mxDOUBLE_CLASS,mxREAL);
    oimy = mxGetPr(plhs[1]);
    
    plhs[3] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
                                   mxGetDimensions(prhs[0]),mxDOUBLE_CLASS,mxREAL);
    oimz = mxGetPr(plhs[3]);
    
    /* Execute the distance transform */
    do_it(iima,dim,oimd,oimx,oimy,oimz);
}
