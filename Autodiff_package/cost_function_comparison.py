import matplotlib.pyplot as plt
import numpy as np
import torch
import time, timeit


def mask_term_gradient(mask_M, mask_omega):
    """
    The mask term of the symbolic gradient expressions 

    Parameters: 
    -----------
    mask_M: np.ndarray
        Binary mask found from optimising a0 to exclude points that are not a part of the sample. 
    mask_omega: np.ndarray
        Binary mask for all valid data points. It is zero for bad detector angular sectors etc.
    """
    return 4* mask_M * mask_omega

def intensity_term_gradient(estimated_I, measured_I, transparency):
    """
    The intensity term of the symbolic gradient expressions 

    Parameters: 
    -----------
    estimated_I: np.ndarray
        The calculated intensity from the current combination of spherical harmonics at the current scanning point
    measured_I: np.ndarray
        The analytical measured intensity
    transparancy: np.ndarray
        The measured intensity needs to account for the transparancy of the sample.
    """
    return ( estimated_I**0.5 - (measured_I/transparency)**0.5) / estimated_I**0.5

def sum_of_spherical_harmonics(a_coeffs, Y_lms):
    """
    The sum of spherical harmonics used in the symbolic gradient expressions

    Parameters: 
    -----------
    a_coeffs: np.ndarray
        All current spherical harmonics coefficients
    Y_lms: np.ndarray
        All current spherical harmonics functions

    """
    return np.sum(a_coeffs * Y_lms, axis = (3,4) ) # Axis order y, x, z, l, m

def epsilon_coeff(mask_M, mask_omega, estimated_I, measured_I, transparency, Y_lm, a_coeffs, Y_lms ):
    """
    Symbolic expression of the gradient to the original cost function with respect to a spherical harmonic coefficient

    Parameters: 
    -----------
    mask_M: np.ndarray
        Binary mask found from optimising a0 to exclude points that are not a part of the sample. 
    mask_omega: np.ndarray
        Binary mask for all valid data points. It is zero for bad detector angular sectors etc.
    estimated_I: np.ndarray
        The calculated intensity from the current combination of spherical harmonics at the current scanning point
    measured_I: np.ndarray
        The analytical measured intensity
    transparancy: np.ndarray
        The measured intensity needs to account for the transparancy of the sample.
    Y_lm: np.ndarray
        The spherical harmonics function belonging to the coefficient that is being optimised
    a_coeffs: np.ndarray
        All current spherical harmonics coefficients
    Y_lms: np.ndarray
        All current spherical harmonics functions
    
    Returns:
        The gradient with respect to the given spherical harmonics coefficient shape y,x,z,l,m
    """

    eps_q = mask_term_gradient(mask_M, mask_omega) * intensity_term_gradient(estimated_I, measured_I, transparency)[:, :, :] # n, y, x, (q,m)
    eps_q = np.sum(eps_q, axis = (0,3)) # Sum over n and (q,phi) -> y,x
    eps_q = eps_q * Y_lm[:,:]
    eps_q = eps_q * sum_of_spherical_harmonics(a_coeffs, Y_lms)
    return eps_q

def torch_estimated_I(a_coeffs, Y_lms)-> torch.tensor:
    """
    Estimate Intensity of one projection like eq (3) in Small-angle X-ray scattering tensor tomography...
    
    Parameters: 
    -----------
    a_coeffs: torch.tensor
        Spherical Harmonics coefficients of voxels with shape (y,x,z,l,m) 
    Y_lms: torch.tensor
        Spherical Harmonics function of voxels with shape (y,x,z,l,m)
    
    Returns: torch.tensor
        Estimated intensity for a given projection with shape (y,x,8). The 8 is the number of segments of radial integration in q-space. 
    """
    I_n = torch.sum( a_coeffs * Y_lms, dim= [3,4] ) 
    I_n = torch.abs( I_n )**2
    return torch.sum( I_n, dim = 2 ) 


def torch_cost_function( mask_omega, Y_lms, measured_I, transparency, a_coeffs)-> torch.tensor:

    estimated_I = torch_estimated_I(a_coeffs, Y_lms)
    eps_q = 2* torch.sum( mask_omega * ( torch.sqrt( estimated_I ) - torch.sqrt( measured_I / transparency) )**2 )
    return eps_q










# Turned quite complicated quite fast. Not necessarly necessary...
# def epsilon_theta_op(mask_M, mask_omega, estimated_I, measured_I, transparency, Y_lm, a_coeffs, Y_lms):
#     """
#     Symbolic expression of the gradient to the original cost function with respect to a the preffered orientation (theta)

#     Parameters: 
#     -----------
#     mask_M: np.ndarray
#         Binary mask found from optimising a0 to exclude points that are not a part of the sample. 
#     mask_omega: np.ndarray
#         Binary mask for all valid data points. It is zero for bad detector angular sectors etc.
#     estimated_I: np.ndarray
#         The calculated intensity from the current combination of spherical harmonics at the current scanning point
#     measured_I: np.ndarray
#         The analytical measured intensity
#     transparancy: np.ndarray
#         The measured intensity needs to account for the transparancy of the sample.
#     Y_lm: np.ndarray
#         The spherical harmonics function belonging to the coefficient that is being optimised
#     a_coeffs: np.ndarray
#         All current spherical harmonics coefficients
#     Y_lms: np.ndarray
#         All current spherical harmonics functions
    
#     Returns:
#         The gradient with respect to the theta_op
#     """
#     eps_q = mask_term_gradient(mask_M, mask_omega)
#     eps_q *= intensity_term_gradient(estimated_I, measured_I, transparency)
#     eps_q = np.sum(eps_q, axis = (0,4)) # Sum over n and phi
#     eps_q *= # Enter function for eq 12 etc.
#     eps_q *= sum_of_spherical_harmonics(a_coeffs, Y_lms)
#     return eps_q


