
import numpy as np
import torch
from initial_comparisons import AD_torch_general, numerical_diff_general
from cost_function_comparison import torch_cost_function, epsilon_coeff, torch_estimated_I
import time



mask_omega = torch.ones( size = (1,1,1) ) # This is 3D in the order: y, x, z. Describes voxels

a_coeffs = torch.ones(size = (1,1,1,1,1) ) *69  # y, x, z, l, m
Y_lms = torch.ones(size = (1,1,1,1,1) )  *169   # y, x, z, l, m ,(8) # Where the last eight represents the shape of the spherical harmoncis 

# Hence, estimated intensity, sums over z, l, m. 

estimated_Is = torch.ones( size = (1,1,1,1 ) ) # n, y, x, (q,phi)=8 ; Goes from voxels to projections. Project a 2D scattering image. As determined in cell above. However radial integration into 8 intensity points. 
estimated_Is = torch_estimated_I(a_coeffs, Y_lms)
measured_Is  = torch.ones( size = (1,1,1,1 ) ) *19 # i.e. projections. 
transparency = torch.ones( size = (1,1,1) ) * 0.5


# ISSUE:
# Going from each voxel to all estimated intensity is not so straight forward, as it is conducted by utilising a range of different functions, coordinate systems, and rotations.
# Moreover, it is somewhat unclear to me how one coefficient and one spherical harmonics function shall represent 8 intensity segments. 
# Seems like there is something about retrieving Y_lm coeff as vectors and matrix multiply with powers of cosine.

# Testing purposes, assume one intensity

print(a_coeffs)

times = 10000

tic = time.time()
for i in range(times):
    grad = AD_torch_general(torch_cost_function, mask_omega, Y_lms, measured_Is, transparency , a_coeffs= a_coeffs) 
tac = time.time()


print(f"Torch autograd: {grad} in {tac-tic} s \n")

tic = time.time()
for i in range(times):
    num_grad = numerical_diff_general(torch_cost_function, mask_omega, Y_lms, measured_Is, transparency , a_coeffs= a_coeffs )
tac = time.time()
print(f"Numerical diff: {num_grad} in {tac-tic} s \n")

tic = time.time()
for i in range(times):
    sym_grad = epsilon_coeff(mask_omega.numpy(), mask_omega.numpy(), estimated_Is.numpy(), measured_Is.numpy(), transparency.numpy(), Y_lms.numpy(), a_coeffs.numpy(), Y_lms.numpy() )
tac = time.time()
print(f"Symbolic gradient: {sym_grad} in {tac-tic} s \n")


"""
Next question is how to introduce 8 intensities. 
The answer, I believe, is that each intensity In(q,phi) is describing one of the 8 segments for a given scanning point(x,y)
We can only find the superposition of spherical harmonics that minimize the error across all projections. 
A given superposition of spherical harmonics will give a certain intensity value across a scanning point with a certain orientation.
The l and m order of spherical harmonics has to mean something for the 8 segments... How to express that? In matlab, they represent a vector...
"""




