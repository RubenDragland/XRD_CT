/*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2016 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|      Authors: Marianne Liebi (marianne.liebi@psi.ch)                  |
%|               Manuel Guizar-Sicairos (manuel.guizar-sicairos@psi.ch)  |
%|               Ivan Usov (ivan.usov@psi.ch)                            |
%|               Oliver Bunk (oliver.bunk@psi.ch)                        |
%|                                                                       |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form and additionally you should cite:
%     M. Liebi, M. Georgiadis, A. Menzel, P. Schneider, J. Kohlbrecher, 
%     O. Bunk, and M. Guizar-Sicairos, “Nanostructure surveys of 
%     macroscopic specimens by small-angle scattering tensor tomography,”
%     Nature 527, 349-352 (2015).   (doi:10.1038/nature16056)
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

#include "mex.h"

void mexFunction(int n_out, mxArray *p_out[], int n_in, const mxArray *p_in[])
{
    double *tomo_obj; // inputs
    double *Ax, *Ay, *Tx, *Ty;
    double temp1, temp2, temp3, temp4;
    size_t nElements, shift_in, shift_out;
    
    double *proj_out; // output
    mwSize *size_out; // size of the output
    size_t nRows, nCols, nPages;
    size_t page_in, page_out; // number of elements per page
    
    size_t i, j, ind; // utility indexes
    
    // read input data
    size_out = (mwSize *) mxGetData(p_in[0]);
    nRows = size_out[0];
    nCols = size_out[1];
    nPages = size_out[2];
    page_out = nRows*nCols;
    
    tomo_obj = mxGetPr(p_in[1]);
    
    page_in = mxGetM(p_in[2])*mxGetN(p_in[2]);
    nElements = mxGetNumberOfElements(p_in[2]);
    Ax = mxGetPr(p_in[2]);
    Ay = mxGetPr(p_in[3]);
    Tx = mxGetPr(p_in[4]);
    Ty = mxGetPr(p_in[5]);
    
    // allocate memory for the output data, initialized to zero
    p_out[0] = mxCreateNumericArray(3, size_out, mxDOUBLE_CLASS, mxREAL);
    proj_out = mxGetPr(p_out[0]);
    
    // sum up the projected data
    for (i = 0; i < nElements; i++) {
        if ((Ax[i] > 0) && (Ax[i] < nCols) && (Ay[i] > 0) && (Ay[i] < nRows)) {
            ind = (size_t) Ay[i] + nRows*Ax[i];
            
            temp1 = Tx[i]*Ty[i];
            temp2 = Tx[i]*(1-Ty[i]);
            temp3 = (1-Tx[i])*Ty[i];
            temp4 = (1-Tx[i])*(1-Ty[i]);
            
            for (j = 0; j < nPages; j++) {
                shift_in = i + j*page_in;
                shift_out = ind + j*page_out;
                
                proj_out[shift_out] += tomo_obj[shift_in] * temp1;
                proj_out[shift_out-1] += tomo_obj[shift_in] * temp2;
                proj_out[shift_out-nRows] += tomo_obj[shift_in] * temp3;
                proj_out[shift_out-nRows-1] += tomo_obj[shift_in] * temp4;
            }
        }
    }
}

