#include "mex.h"
#include <math.h>


void add_3d(double *mat, double *dims, long c1, long c2, long c3, double toadd){
    if ((c1>0)&&(c1<=dims[0])&&(c2>0)&&(c2<=dims[1])&&(c3>0)&&(c3<=dims[2]))
        mat[(c1-1)+(c2-1)*(int)dims[0]+(c3-1)*(int)dims[0]*(int)dims[1]]+=toadd;
    return;
}

double get_2d(double *mat, double *dims, long c1, long c2){
    if ((c1>0)&&(c1<=dims[0])&&(c2>0)&&(c2<=dims[1]))
        return mat[(c1-1)+(c2-1)*(int)dims[0]];
    else
        return 0;
}

void mexFunction(int n_out, mxArray *p_out[], int n_in, const mxArray *p_in[]) 
/*(tomo,[startx0,starty0,startz],[stepx,stepy,stepz],[steppixx,steppixy,steppixz],BTranspose[6,8],diffmap[8,106,70])*/
{
    double xp,yp,r;
    int i,j,k,t,xp0,yp0;
    double *proj_out;
    double *tomo, *modelsize, *projsize, *Xp, *Yp, *diffmap;
    double diff;
    
    tomo = mxGetPr(p_in[0]);
    
    modelsize = mxGetPr(p_in[1]);
    r=(double)modelsize[1]/2;
    
    projsize = mxGetPr(p_in[2]);
    
    Xp = mxGetPr(p_in[3]);
    
    Yp = mxGetPr(p_in[4]);
         
    diffmap = mxGetPr(p_in[5]);
    
    mwSize outdims[3];
    
    double *tomonew;
    /*outdims = mxGetDimensions(p_in[0]);*/
    outdims[0]=modelsize[0]; outdims[1]=modelsize[1]; outdims[2]=modelsize[2];
    p_out[0] = mxCreateNumericArray(3, outdims, mxDOUBLE_CLASS, mxREAL);
    tomonew = mxGetPr(p_out[0]);
    
    for (i=0;(i<modelsize[0]*modelsize[1]*modelsize[2]);i++){
        tomonew[i] = tomo[i];
    }
    for (i=1;(i<=modelsize[0]);i++) {
        for (j=1;(j<=modelsize[1]);j++) {
            for (k=1;(k<=modelsize[2]);k++) {
                if ((((double)j-r)*((double)j-r)+((double)k-r)*((double)k-r))<r*r) {
                    xp=Xp[i+j*(int)modelsize[0]+k*(int)modelsize[0]*(int)modelsize[1]];
                    yp=Yp[i+j*(int)modelsize[0]+k*(int)modelsize[0]*(int)modelsize[1]];
                    xp0=(int)floor(xp);
                    yp0=(int)floor(yp);
                    diff= get_2d(diffmap,projsize,xp0,yp0)*((double)xp0+1-xp)*((double)yp0+1-yp) + \
                            get_2d(diffmap,projsize,xp0+1,yp0)*(xp-(double)xp0)*((double)yp0+1-yp) + \
                            get_2d(diffmap,projsize,xp0,yp0+1)*((double)xp0+1-xp)*(yp-(double)yp0) + \
                            get_2d(diffmap,projsize,xp0+1,yp0+1)*(xp-(double)xp0)*(yp-(double)yp0);
                
                    add_3d(tomonew,modelsize,i,j,k,diff);
            
                }
            }
        }
    }
}
