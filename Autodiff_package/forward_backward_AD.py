import torch
import numpy as np






def forward_backward(cost_function, a_temp: torch.tensor, Ylm, voxel_dims, numOfsegments, data, current_projection, Rot_exp_now, p, find_grad, X, Y, Z, numOfpixels):
    """
    Python translation of forward pass in reconstruction algorithm.
    Calculates error in reciprocal intensity given the estimated SH functions.
    
    """

    ny, nx, nz = voxel_dims

    if find_grad:
        a_temp.requires_grad_(True) # Assume transormed to torch. May have to adapt AD_torch_general from initial_comparisons.


    error_norm, aux_diff_poisson, proj_out_all  = cost_function( a_temp, Ylm, voxel_dims, numOfsegments, data, current_projection, Rot_exp_now, p, find_grad, X, Y, Z, numOfpixels ) 

    if find_grad:
        error_norm.backward() # No need for args?
        ad_grad = a_temp.grad
    else:
        ad_grad = None # Something like this. Or separate into forward function and forward_backward function. Or only do forward_backward in python.

    return error_norm, ad_grad, aux_diff_poisson, proj_out_all


def format_matlab_input():
    """
    Deals with matlab engine things.
    After this, assume everything is in torch.
    """
    return


def SAXS_AD_cost_function(a_temp, Ylm, ny, nx, nz, numOfsegments, data, current_projection, Rot_exp_now, p, find_grad, X, Y, Z, numOfpixels):
    """
    Cost function of coefficients AD-SAXS
    """

    sumlm_alm_Ylm = page_multiply(a_temp, Ylm)

    data_synt_vol = torch.permute( torch.abs(sumlm_alm_Ylm**2), (3,2,1) )
    data_synt_vol = torch.reshape( data_synt_vol, (ny,nx,nz, numOfsegments) )

    xout = np.arange(1, data.size(1) ) - torch.ceil( data.size(1)/2 ) + current_projection["dx"] # NB indices
    yout = np.arange(1, data.size(0) ) - torch.ceil( data.size(0)/2 ) + current_projection["dy"]

    proj_out_all = arb_projection(data_synt_vol, X, Y, Z, Rot_exp_now, p, xout, yout)

    aux_diff_poisson = ( torch.sqrt(proj_out_all) - torch.sqrt(data) ) * current_projection["window_mask"]

    error_norm = 2 * torch.sum( aux_diff_poisson**2, dim = (2,1,0)) # Check dims

    return error_norm, aux_diff_poisson, proj_out_all




def page_multiply(A: torch.tensor, B: torch.tensor): # Check if numba compatible

    assert (A.size[-1] == B.size[-1] )

    for i in range( A.size[-1] ):

        A[:,:,i] = torch.matmul(A[:,:,i],B[:,:,i]) # Check if works compared to creating new Tensor.

    return A 

def arb_projection(tomo_obj_all, X, Y, Z, R, p, xout, yout): # Assume numpy
    """
    Estimate projection based on SH. Ignore most options for now.
    """
    try:
        method = p["method"] # Retrieve method and mode
        mode = p["mode"]
    except:
        method = "nearest"
        mode = 1
    finally:
    

        Xp = R[0,0]*X + R[0,1]*Y + R[0,2]*Z  # Rot matrix is 3x3. X, Y, Z are meshgrid holding indices. NB! 0 indexing python vs 1 indexing matlab.
        Yp = R[1,0]*X + R[1,1]*Y + R[1,2]*Z 

        # Volume upsampling ignored

        min_xout = np.min(xout)
        min_yout = np.min(yout)

        Ax = np.round(Xp-min_xout) # Drop +1 due to 0-indexing
        Ay = np.round(Yp-min_yout) 

        proj_out_all = torch.zeros( len(yout), len(xout), tomo_obj_all)



    return


def nearest_interpolation():

    return


def bilinear_interpolation():

    return
