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

// Wrapper function
void mexFunction(int n_out, mxArray *p_out[], int n_in, const mxArray *p_in[])
{
    double *in; // input image
    const mwSize *size_in; // size of input image
    
    double *out; // interpolated image
    mwSize *size_out; // size of interpolated image
    size_t nRows, nCols, nPages, nSeries; // dimentions of interpolated image
    
    mwSize ndim; // number of dimensions for both images (3 or 4)
    
    size_t i, j, k, l; // utility indexes
    size_t page, series; // utility variables
    
    // check number of inputs
    if (n_in != 1) {
        mexErrMsgIdAndTxt("linupsamp3:invalidInputsNum",
                "Function takes 1 input arrays");
    }
    
    // check type of inputs
    if (!mxIsDouble(p_in[0])) {
        mexErrMsgIdAndTxt("linupsamp3:invalidInputType",
                "Function takes input array of type double");
    }
    
    // check number of dimensions
    ndim = mxGetNumberOfDimensions(p_in[0]);
    if (ndim < 3 || 4 < ndim) {
        mexErrMsgIdAndTxt("linupsamp3:invalidInputDimensions",
                "Function takes input 3D or 4D arrays");
    }
    
    // read input image
    in = mxGetPr(p_in[0]);
    size_in = mxGetDimensions(p_in[0]);
    
    // assign dimensions of the output array
    size_out = (mwSize *) mxMalloc(ndim * sizeof(mwSize));
    for (i = 0; i < 3; i++) {
        size_out[i] = 2*size_in[i]-1;
    }
    
    if (ndim > 3) {
        size_out[3] = size_in[3];
    }
    
    nRows = size_out[0];
    nCols = size_out[1];
    nPages = size_out[2];
    nSeries = (ndim > 3) ? size_out[3] : 1;
    
    // create output image array
    p_out[0] = mxCreateNumericArray(ndim, size_out, mxDOUBLE_CLASS, mxREAL);
    out = mxGetPr(p_out[0]);
    
    // free memory allocated by mxMalloc
    mxFree(size_out);
    
    // interpolate
    page = nRows*nCols;
    series = page*nPages;
    
    // copy vertex points from the initial image
    size_t series_shift, page_shift, col_shift;
    size_t val = 0;
    for (i = 0; i < nSeries; i++) {
        series_shift = i*series;
        
        for (j = 0; j < nPages; j+=2) {
            page_shift = series_shift + j*page;
            
            for (k = 0; k < nCols; k+=2) {
                col_shift = page_shift + k*nRows;
                
                for (l = 0; l < nRows; l+=2) {
                    out[col_shift+l] = in[val++];
                }
            }
        }
    }
    
    // values within rows
    for (i = 0; i < nSeries; i++) {
        series_shift = i*series;
        
        for (j = 0; j < nPages; j+=2) {
            page_shift = series_shift + j*page;
            
            for (k = 0; k < nCols; k+=2) {
                col_shift = page_shift + k*nRows;
                
                for (l = 1; l < nRows; l+=2) {
                    val = col_shift+l;
                    out[val] = (out[val-1] + out[val+1]) / 2.0;
                }
            }
        }
    }
    
    // values within planes
    for (i = 0; i < nSeries; i++) {
        series_shift = i*series;
        
        for (j = 0; j < nPages; j+=2) {
            page_shift = series_shift + j*page;
            
            for (k = 1; k < nCols; k+=2) {
                col_shift = page_shift + k*nRows;
                
                for (l = 0; l < nRows; l++) {
                    val = col_shift+l;
                    out[val] = (out[val-nRows] + out[val+nRows]) / 2.0;
                }
            }
        }
    }
    
    // values within volumes
    for (i = 0; i < nSeries; i++) {
        series_shift = i*series;
        
        for (j = 1; j < nPages; j+=2) {
            page_shift = series_shift + j*page;
            
            for (k = 0; k < nCols; k++) {
                col_shift = page_shift + k*nRows;
                
                for (l = 0; l < nRows; l++) {
                    val = col_shift+l;
                    out[val] = (out[val-page] + out[val+page]) / 2.0;
                }
            }
        }
    }
}

