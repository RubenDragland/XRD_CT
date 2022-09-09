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
#include "blas.h"

// Wrapper function
void mexFunction(int n_out, mxArray *p_out[], int n_in, const mxArray *p_in[])
{
    double *A, *B, *C; // pointers to input and output real parts
    double *Ai, *Bi, *Ci; // pointers to input and output imaginary part
    bool isComplexA, isComplexB; // complex input flags
    mwSize ndimA, ndimB, ndimC; // number of dimensions
    const mwSize *sizeA, *sizeB; // sizes of input matrices
    size_t rowA, colA, rowB, colB; // number of rows and columns
    mwSize *sizeC; // size of output matrix
    
    bool bsx = false; // flag for singleton expantion
    mwSize *sizeA_pad, *sizeB_pad; // padded sizes of matrices after singleton expansion
    int *pageIndA, *pageIndB; // page indexes of each multiplication
    
    size_t pageA, pageB, pageC; // page size (page = row * col)
    double *A_page, *B_page, *C_page; // pointers to pages of real parts
    double *Ai_page, *Bi_page, *Ci_page; // pointers to pages of imaginary parts
    
    mwSize *idx, indA, indB; // utility indexes
    int i, j;
    
    // check number of inputs
    if (n_in != 2) {
        mexErrMsgIdAndTxt("bsxpagemult:invalidNumInputs",
                "Function takes 2 input arrays");
    }
    
    // check type of inputs
    if (!mxIsDouble(p_in[0]) || !mxIsDouble(p_in[1])) {
        mexErrMsgIdAndTxt("bsxpagemult:invalidTypeInputs",
                "Function takes input arrays of type double");
    }
    
    isComplexA = mxIsComplex(p_in[0]);
    isComplexB = mxIsComplex(p_in[1]);
    
    if (isComplexA && isComplexB) {
        mexErrMsgIdAndTxt("bsxpagemult:invalidTypeInputs",
                "Both input arrays cannot be complex");
    }
    
    // read input matrix A
    A = mxGetPr(p_in[0]);
    if (isComplexA) {
        Ai = mxGetPi(p_in[0]);
    }
    ndimA = mxGetNumberOfDimensions(p_in[0]);
    sizeA = mxGetDimensions(p_in[0]);
    rowA = sizeA[0];
    colA = sizeA[1];
    
    // make sure that any dimension of A is not equal to 0
    for (i = 0; i < ndimA; i++) {
        if (sizeA[i] == 0) {
            mexErrMsgIdAndTxt("bsxpagemult:unsupportedZeroDimension",
                    "Matrix A has unsupported zero dimension(s)");
        }
    }
    
    // read input matrix B
    B = mxGetPr(p_in[1]);
    if (isComplexB) {
        Bi = mxGetPi(p_in[1]);
    }
    ndimB = mxGetNumberOfDimensions(p_in[1]);
    sizeB = mxGetDimensions(p_in[1]);
    rowB = sizeB[0];
    colB = sizeB[1];
    
    // make sure that any dimension of B is not equal to 0
    for (i = 0; i < ndimB; i++) {
        if (sizeB[i] == 0) {
            mexErrMsgIdAndTxt("bsxpagemult:unsupportedZeroDimension",
                    "Matrix B has unsupported zero dimension(s)");
        }
    }
    
    if (colA != rowB) {
        mexErrMsgIdAndTxt("bsxpagemult:invalidSizeInputs",
                "Number of columns in A must be equal to number of rows in B");
    }
    
    // process outputs
    ndimC = (ndimA > ndimB) ? ndimA : ndimB;
    sizeC = (mwSize *) mxMalloc(ndimC * sizeof(mwSize));
    
    sizeC[0] = rowA; // rowC is equal to rowA
    sizeC[1] = colB; // colC is equal to colB
    
    sizeA_pad = (mwSize *) mxMalloc(ndimC * sizeof(mwSize));
    sizeB_pad = (mwSize *) mxMalloc(ndimC * sizeof(mwSize));
    
    // Pad input array sizes and calculate the total number of multiplications
    int nMult = 1;
    for (i = 0; i < ndimC; i++) {
        sizeA_pad[i] = (i < ndimA) ? sizeA[i] : 1;
        sizeB_pad[i] = (i < ndimB) ? sizeB[i] : 1;
        
        if (i > 1) {
            if (sizeA_pad[i] != 1 && sizeB_pad[i] != 1 && sizeA_pad[i] != sizeB_pad[i]) {
                mexErrMsgIdAndTxt("bsxpagemult:invalidSizeInputs",
                        "Non-singleton dimensions (above 2) of inputs must match each other.");
            }
            
            if (sizeA_pad[i] != sizeB_pad[i]) {
                bsx = true;
            }
            
            sizeC[i] = (sizeA_pad[i] > sizeB_pad[i]) ? sizeA_pad[i] : sizeB_pad[i];
            nMult *= sizeC[i];
        }
    }
    
    // create output matrix C
    if (isComplexA || isComplexB) {
        p_out[0] = mxCreateUninitNumericArray(ndimC, sizeC, mxDOUBLE_CLASS, mxCOMPLEX);
        C = mxGetPr(p_out[0]);
        Ci = mxGetPi(p_out[0]);
    }
    else {
        p_out[0] = mxCreateUninitNumericArray(ndimC, sizeC, mxDOUBLE_CLASS, mxREAL);
        C = mxGetPr(p_out[0]);
    }
    
    // compute page indexes for each multiplication
    if (bsx) {
        idx = (mwSize *) mxMalloc(ndimC * sizeof(mwSize));
        for (i = 0; i < ndimC; i++) {
            idx[i] = 0;
        }
        
        pageIndA = (int *) mxMalloc(nMult * sizeof(int));
        pageIndB = (int *) mxMalloc(nMult * sizeof(int));
        pageIndA[0] = pageIndB[0] = 0;
        
        for (i = 1; i < nMult; i++) {
            // idx = ind2sub(size(C), i) - 1
            idx[2]++;
            for (j = 2; j < ndimC; j++) {
                if (idx[j] > sizeC[j]-1) {
                    idx[j] = 0;
                    idx[j+1]++;
                }
            }
            
            // indA = sub2ind(sizeA_pad, idx + 1) - 1
            // indB = sub2ind(sizeB_pad, idx + 1) - 1
            indA = indB = 0;
            for (j = ndimC-1; j > 1; j--) {
                if (sizeA_pad[j] > 1) {
                    indA = indA*sizeA_pad[j] + idx[j];
                }
                
                if (sizeB_pad[j] > 1) {
                    indB = indB*sizeB_pad[j] + idx[j];
                }
            }
            
            // assign
            pageIndA[i] = indA;
            pageIndB[i] = indB;
        }
    }
    
    pageA = rowA*colA;
    pageB = rowB*colB;
    pageC = rowA*colB;
    
    // loop over page matrix multiplications
    double one = 1.0, zero = 0.0;
    char *ch = "N";
    
    if (isComplexA) {
        for (i = 0; i < nMult; i++) {
            C_page = C + i*pageC;
            Ci_page = Ci + i*pageC;
            
            if (bsx) {
                A_page = A + pageIndA[i]*pageA;
                Ai_page = Ai + pageIndA[i]*pageA;
                B_page = B + pageIndB[i]*pageB;
            }
            else {
                A_page = A + i*pageA;
                Ai_page = Ai + i*pageA;
                B_page = B + i*pageB;
            }
            
            dgemm(ch, ch, &rowA, &colB, &colA,
                    &one, A_page, &rowA, B_page, &rowB,
                    &zero, C_page, &rowA);
            
            dgemm(ch, ch, &rowA, &colB, &colA,
                    &one, Ai_page, &rowA, B_page, &rowB,
                    &zero, Ci_page, &rowA);
        }
    }
    else if (isComplexB) {
        for (i = 0; i < nMult; i++) {
            C_page = C + i*pageC;
            Ci_page = Ci + i*pageC;
            
            if (bsx) {
                A_page = A + pageIndA[i]*pageA;
                B_page = B + pageIndB[i]*pageB;
                Bi_page = Bi + pageIndB[i]*pageB;
            }
            else {
                A_page = A + i*pageA;
                B_page = B + i*pageB;
                Bi_page = Bi + i*pageB;
            }
            
            dgemm(ch, ch, &rowA, &colB, &colA,
                    &one, A_page, &rowA, B_page, &rowB,
                    &zero, C_page, &rowA);
            
            dgemm(ch, ch, &rowA, &colB, &colA,
                    &one, A_page, &rowA, Bi_page, &rowB,
                    &zero, Ci_page, &rowA);
        }
    }
    else {
        for (i = 0; i < nMult; i++) {
            C_page = C + i*pageC;
            
            if (bsx) {
                A_page = A + pageIndA[i]*pageA;
                B_page = B + pageIndB[i]*pageB;
            }
            else {
                A_page = A + i*pageA;
                B_page = B + i*pageB;
            }
            
            dgemm(ch, ch, &rowA, &colB, &colA,
                    &one, A_page, &rowA, B_page, &rowB,
                    &zero, C_page, &rowA);
        }
    }
    
    // free memory allocated by mxMalloc
    mxFree(sizeC);
    mxFree(sizeA_pad);
    mxFree(sizeB_pad);
    if (bsx) {
        mxFree(idx);
        mxFree(pageIndA);
        mxFree(pageIndB);
    }
}

