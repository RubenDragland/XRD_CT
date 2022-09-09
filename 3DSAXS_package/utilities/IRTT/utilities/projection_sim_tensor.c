// projection_sim_tensor(tomotensor,modelsize,projsize,Xp,Yp,B)
// Mex function, calling it without pipeline might cause matlab to crash. Use proj_simul_tensor() instead.
// To compile with openMP parallel processing, do " mex 'CFLAGS="\$CFLAGS -std=c99 -fopenmp"' LDFLAGS="\$LDFLAGS -fopenmp" *.c ".
// " mex *.c " to compile without parallel.

// Copyright 2017 Zirui Gao
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "mex.h"
#include <math.h>

void add_3d(double *mat, long *dims, long c1, long c2, long c3, double toadd){
    if ((c1>0)&&(c1<=dims[0])&&(c2>0)&&(c2<=dims[1])&&(c3>0)&&(c3<=dims[2]))
        mat[(c1-1)+(c2-1)*dims[0]+(c3-1)*dims[0]*dims[1]]+=toadd;
    return;
}

void mexFunction(int n_out, mxArray *p_out[], int n_in, const mxArray *p_in[]) 
/*(tomo,[startx0,starty0,startz],[stepx,stepy,stepz],[steppixx,steppixy,steppixz],Binv[6,(long)projsize[2]])*/
{
    double xp,yp,r,temp;
    long i,j,k,l,t,xp0,yp0,size_loop,loop;
    double *proj_out;
    double *tomo, *modelsize, *projsize, *Xp, *Yp, *B;
    
    tomo = mxGetPr(p_in[0]);
    
    modelsize = mxGetPr(p_in[1]);
    long msize[4];
    for (i=0;i<=3;i++) msize[i]=(long)modelsize[i];
    r=modelsize[1]/2;
    
    projsize = mxGetPr(p_in[2]);
    long psize[3];
    for (i=0;i<=2;i++) psize[i]=(long)projsize[i];
    
    Xp = mxGetPr(p_in[3]);
    
    Yp = mxGetPr(p_in[4]);
     
    B = mxGetPr(p_in[5]);
    
    mwSize outdims[3];
    outdims[0] = psize[0]; outdims[1] = psize[1]; outdims[2] = psize[2];
    p_out[0] = mxCreateNumericArray(3, outdims, mxDOUBLE_CLASS, mxREAL);
    proj_out = mxGetPr(p_out[0]);
    
    size_loop=msize[0]*msize[1]; // For openmp parallel processing
        
    #pragma omp parallel private(i,j,k,l,t,xp,yp,xp0,yp0,temp)
    #pragma omp for
    for (loop=0; loop< size_loop; loop++) {
        i=loop/msize[1];
        j=loop-i*msize[1];
        for (k=0;k<=(msize[2]-1);k++){
            if ((((double)j-r)*((double)j-r)+((double)k-r)*((double)k-r))<r*r) {
                xp=Xp[i+j*msize[0]+k*msize[0]*msize[1]];
                yp=Yp[i+j*msize[0]+k*msize[0]*msize[1]];
                xp0=(long)floor(xp);
                yp0=(long)floor(yp);
                if ((xp>=1) && (xp<(psize[1])) && (yp>=1) && (yp<(psize[2]))) {
                    for (l=1;l<=psize[0];l++){
                        temp=0;
                        for (t=0;t<=(msize[3]-1);t++){
                            temp+=tomo[i+j*msize[0]+k*msize[0]*msize[1]+t*msize[0]*msize[1]*msize[2]]*B[(l-1)+t*psize[0]];
                        }
#pragma omp atomic
                        proj_out[(l-1)+(xp0-1)*psize[0]+(yp0-1)*psize[0]*psize[1]] += temp * ((double)xp0+1-xp) * ((double)yp0+1-yp);
#pragma omp atomic
                        proj_out[(l-1)+xp0*psize[0]+(yp0-1)*psize[0]*psize[1]] += temp * (xp-(double)xp0) * ((double)yp0+1-yp);
#pragma omp atomic
                        proj_out[(l-1)+(xp0-1)*psize[0]+yp0*psize[0]*psize[1]] += temp * ((double)xp0+1-xp) * (yp-(double)yp0);
#pragma omp atomic
                        proj_out[(l-1)+xp0*psize[0]+yp0*psize[0]*psize[1]] += temp * (xp-(double)xp0) * (yp-(double)yp0);
                    }
                }
            }
        }
    }
    
}
