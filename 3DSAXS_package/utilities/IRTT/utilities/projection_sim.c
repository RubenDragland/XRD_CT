#include "mex.h"
#include <math.h>

void add_2d(double *mat, double *dims, long c1, long c2, double toadd){
    if ((c1>0)&&(c1<=dims[0])&&(c2>0)&&(c2<=dims[1]))
        mat[(c1-1)+(c2-1)*(int)dims[0]]+=toadd;
    return;
}

void mexFunction(int n_out, mxArray *p_out[], int n_in, const mxArray *p_in[])
{
    double xp,yp;
    int i,j,k,t,xp0,yp0;
    double *proj_out;
    double *tomo, *modelsize, *projsize, *Xp, *Yp;
    
    tomo = mxGetPr(p_in[0]);
    
    modelsize = mxGetPr(p_in[1]);
    
    projsize = mxGetPr(p_in[2]);
    
    Xp = mxGetPr(p_in[3]);
    
    Yp = mxGetPr(p_in[4]);
        
    p_out[0] = mxCreateDoubleMatrix(projsize[0], projsize[1], mxREAL);
    proj_out = mxGetPr(p_out[0]);
    
    for (i=0;i<=(modelsize[0]-1);i++) {
        for (j=0;j<=(modelsize[1]-1);j++) {
            for (k=0;k<=(modelsize[2]-1);k++){
                xp=Xp[i+j*(int)modelsize[0]+k*(int)modelsize[0]*(int)modelsize[1]];
                yp=Yp[i+j*(int)modelsize[0]+k*(int)modelsize[0]*(int)modelsize[1]];
                //if ((xp>=1) && (xp<=(projsize[0])) && (yp>=1) && (yp<=(projsize[1]))) {
                    xp0=(int)floor(xp);
                    yp0=(int)floor(yp);
                    add_2d(proj_out,projsize,xp0,yp0,tomo[i+j*(int)modelsize[0]+k*(int)modelsize[0]*(int)modelsize[1]] * ((double)xp0+1-xp) * ((double)yp0+1-yp));
                    add_2d(proj_out,projsize,xp0+1,yp0,tomo[i+j*(int)modelsize[0]+k*(int)modelsize[0]*(int)modelsize[1]] * (xp-(double)xp0) * ((double)yp0+1-yp));
                    add_2d(proj_out,projsize,xp0,yp0+1,tomo[i+j*(int)modelsize[0]+k*(int)modelsize[0]*(int)modelsize[1]] * ((double)xp0+1-xp) * (yp-(double)yp0));
                    add_2d(proj_out,projsize,xp0+1,yp0+1,tomo[i+j*(int)modelsize[0]+k*(int)modelsize[0]*(int)modelsize[1]] * (xp-(double)xp0) * (yp-(double)yp0));
                //}
            }
        }
    }
}
