import matlab.engine
import numpy as np
import sys
import os
import torch


def exp_sin_squared(A, B, SIN_SQUARE_THETA):
    return A**2 * torch.exp(-B * SIN_SQUARE_THETA)


def int_cost_function(A, I_n, n):
    I_n_hat = torch.sum(A, axis=n)  # Sum over the axis given by n. x, y or z.
    return torch.sum((I_n_hat - I_n) ** 2)


def lab_rotation_given_axis_n(n):
    # Returns alpha, beta
    if n == 0:
        return 0, np.pi / 2
    elif n == 1:
        return np.pi / 2, 0
    elif n == 2:
        return 0, 0
    else:
        raise ValueError("n must be 0, 1 or 2.")


def derive_SIN_THETA(
    theta_op,
    phi_op,
    n,
    pos_phi,
):
    """
    Derives the SIN_THETA matrix for the given n and phi.
    n is axis number, and must be 0,1,2.
    phi is spherical position coordinates. theta position is evaluated at pi/2.
    """
    # Derive the Theta function
    alpha, beta = lab_rotation_given_axis_n(n)
    alpha = torch.tensor(alpha)
    beta = torch.tensor(beta)
    pos_phi = torch.tensor(pos_phi)

    # Rotation matrices
    COS_THETA = (
        torch.sin(theta_op) * torch.cos(phi_op) * torch.cos(beta)
        + torch.cos(theta_op) * torch.sin(beta)
    ) * torch.cos(pos_phi)
    +(
        torch.sin(theta_op) * torch.cos(phi_op) * torch.sin(alpha) * torch.cos(beta)
        + torch.sin(theta_op) * torch.sin(phi_op) * torch.cos(alpha)
        - torch.cos(theta_op) * torch.cos(beta) * torch.sin(alpha)
    ) * torch.sin(pos_phi)

    SIN_THETA_SQUARED = 1 - COS_THETA**2

    return SIN_THETA_SQUARED


def new_cost_function(A, B, theta, phi, I_n, n, simple=True, pos_phi=0):
    # alpha, beta = lab_rotation_given_axis_n(n)
    if simple:
        SIN_SQUARE_THETA = torch.sin(theta) ** 2
    else:
        SIN_SQUARE_THETA = derive_SIN_THETA(theta, phi, n, pos_phi=pos_phi)

    # Temp rotation # TODO: Get from matlab. Implement properly... Does not necessarily make sense. But will be the same issue in both implementations
    I_n_hat = torch.sum(exp_sin_squared(A, B, SIN_SQUARE_THETA), axis=n)
    eps_q = 2 * torch.sum((I_n_hat**0.5 - I_n**0.5) ** 2)
    return eps_q


def int_model_gradients(A, I_n, n):
    # Only supports square tensors
    grad_A = torch.zeros(A.shape)
    I_n_hat = torch.sum(A, axis=n)
    for i in range(A.shape[-1]):
        grad_A[:, :, i] = 2 * (I_n_hat - I_n)
    return torch.permute(
        grad_A, (2, 0, 1) if n == 0 else (0, 2, 1) if n == 1 else (0, 1, 2)
    )


def BP(A, B, SST, I_n, n):
    return (
        torch.sum(exp_sin_squared(A, B, SST), axis=n) ** 0.5 - I_n**0.5
    ) / torch.sum(exp_sin_squared(A, B, SST), axis=n) ** 0.5


def new_model_grad_A(A, B, theta, phi, I_n, n, simple=True, pos_phi=0):
    grad_A = torch.zeros(A.shape)
    if simple:
        SIN_SQUARE_THETA = torch.sin(theta) ** 2
    else:
        SIN_SQUARE_THETA = derive_SIN_THETA(theta, phi, n, pos_phi=pos_phi)
    BP_n = BP(A, B, SIN_SQUARE_THETA, I_n, n)
    for i in range(A.shape[-1]):
        # This part is still uncertain... Does it make sense?
        grad_A[:, :, i] = (
            BP_n * exp_sin_squared(A, B, SIN_SQUARE_THETA)[:, :, i] / A[:, :, i]
        )
    # More efficient if index params. Note A**2 term included in exp_sin_squared function. Have to divide by A
    return torch.permute(
        4 * grad_A, (2, 0, 1) if n == 0 else (0, 2, 1) if n == 1 else (0, 1, 2)
    )


def new_model_grad_B(A, B, theta, phi, I_n, n, simple=True, pos_phi=0):
    grad_B = torch.zeros(A.shape)
    if simple:
        SIN_SQUARE_THETA = torch.sin(theta) ** 2
    else:
        SIN_SQUARE_THETA = derive_SIN_THETA(theta, phi, n, pos_phi=pos_phi)

    BP_n = BP(A, B, SIN_SQUARE_THETA, I_n, n)
    for i in range(A.shape[-1]):
        grad_B[:, :, i] = (
            BP_n
            * SIN_SQUARE_THETA[:, :, i]
            * exp_sin_squared(A, B, SIN_SQUARE_THETA)[:, :, i]
        )
        # More efficient if index params. Note A**2 term included in exp_sin_squared function.
    return torch.permute(
        -2 * grad_B, (2, 0, 1) if n == 0 else (0, 2, 1) if n == 1 else (0, 1, 2)
    )  # Check this


def new_model_grad_theta(A, B, theta, phi, I_n, n, simple=True, pos_phi=0):
    """
    Implements the gradient of the cost function with respect to theta.
    Wonder if some wrong assumptions are made.
    """
    grad_theta = torch.zeros(A.shape)
    if simple:
        SIN_SQUARE_THETA = torch.sin(theta) ** 2
    else:
        SIN_SQUARE_THETA = derive_SIN_THETA(theta, phi, n, pos_phi=pos_phi)

    BP_n = BP(A, B, SIN_SQUARE_THETA, I_n, n)
    for i in range(A.shape[-1]):
        grad_theta[:, :, i] = (
            BP_n
            * B[:, :, i]
            * SIN_SQUARE_THETA[:, :, i]
            * exp_sin_squared(A, B, SIN_SQUARE_THETA)[:, :, i]
            * dTHETA_dtheta_op(SIN_SQUARE_THETA, theta, phi, n, pos_phi=pos_phi)[
                :, :, i
            ]
        )

    return -4 * grad_theta


def new_model_grad_phi(A, B, theta, phi, I_n, n, simple=True, pos_phi=0):
    """
    Implements the gradient of the cost function with respect to phi.
    Wonder if some wrong assumptions are made.
    """
    grad_phi = torch.zeros(A.shape)
    if simple:
        SIN_SQUARE_THETA = torch.sin(theta) ** 2
    else:
        SIN_SQUARE_THETA = derive_SIN_THETA(theta, phi, n, pos_phi=pos_phi)

    BP_n = BP(A, B, SIN_SQUARE_THETA, I_n, n)
    for i in range(A.shape[-1]):
        grad_phi[:, :, i] = (
            BP_n
            * B[:, :, i]
            * SIN_SQUARE_THETA[:, :, i]
            * exp_sin_squared(A, B, SIN_SQUARE_THETA)[:, :, i]
            * dTHETA_dphi_op(SIN_SQUARE_THETA, theta, phi, n, pos_phi=pos_phi)[:, :, i]
        )

    return -4 * grad_phi


def dTHETA_dtheta_op(SIN_SQUARE_THETA, theta_op, phi_op, n, pos_phi=0):

    alpha, beta = lab_rotation_given_axis_n(n)
    alpha = torch.tensor(alpha)
    beta = torch.tensor(beta)
    pos_phi = torch.tensor(pos_phi)

    # dcos_dtheta = torch.cos(pos_phi) * (
    #     torch.cos(theta_op) * torch.cos(phi_op) * torch.cos(beta)
    #     - torch.sin(theta_op) * torch.sin(beta)
    # ) + torch.sin(pos_phi) * (
    #     torch.cos(theta_op) * torch.cos(phi_op) * torch.sin(alpha) * torch.cos(beta)
    #     + torch.cos(theta_op) * torch.sin(phi_op) * torch.cos(alpha)
    #     + torch.sin(theta_op) * torch.cos(beta) * torch.sin(alpha)
    # )

    # grad = (1 / torch.sqrt(SIN_SQUARE_THETA)) * dcos_dtheta

    COS_THETA = (
        torch.sin(theta_op) * torch.cos(phi_op) * torch.cos(beta)
        + torch.cos(theta_op) * torch.sin(beta)
    ) * torch.cos(pos_phi)

    x_s = torch.cos(beta) * torch.cos(pos_phi) + torch.sin(alpha) * torch.cos(
        beta
    ) * torch.sin(pos_phi)
    y_s = torch.cos(alpha) * torch.sin(pos_phi)
    z_s = torch.sin(beta) * torch.cos(pos_phi) - torch.cos(beta) * torch.sin(
        alpha
    ) * torch.sin(pos_phi)

    grad = (
        (1 / torch.sin(torch.arccos(COS_THETA)))
        * (
            z_s * torch.sin(theta_op)
            - (x_s * torch.cos(phi_op))
            + y_s * torch.sin(phi_op)
        )
        * torch.cos(theta_op)
    )

    return grad


def dTHETA_dphi_op(SIN_SQUARE_THETA, theta_op, phi_op, n, pos_phi=0):
    alpha, beta = lab_rotation_given_axis_n(n)
    alpha = torch.tensor(alpha)
    beta = torch.tensor(beta)
    pos_phi = torch.tensor(pos_phi)

    # dcos_dphi = torch.cos(pos_phi) * (
    #     -torch.sin(theta_op) * torch.sin(phi_op) * torch.cos(beta)
    # ) + torch.sin(pos_phi) * (
    #     -torch.sin(phi_op) * torch.sin(phi_op) * torch.sin(alpha) * torch.cos(beta)
    #     + torch.sin(theta_op) * torch.cos(phi_op) * torch.cos(alpha)
    # )

    # grad = (1 / torch.sqrt(SIN_SQUARE_THETA)) * dcos_dphi

    COS_THETA = (
        torch.sin(theta_op) * torch.cos(phi_op) * torch.cos(beta)
        + torch.cos(theta_op) * torch.sin(beta)
    ) * torch.cos(pos_phi)

    x_s = torch.cos(beta) * torch.cos(pos_phi) + torch.sin(alpha) * torch.cos(
        beta
    ) * torch.sin(pos_phi)
    y_s = torch.cos(alpha) * torch.sin(pos_phi)
    z_s = torch.sin(beta) * torch.cos(pos_phi) - torch.cos(beta) * torch.sin(
        alpha
    ) * torch.sin(pos_phi)

    grad = (
        (1 / torch.sin(torch.arccos(COS_THETA)))
        * torch.sin(theta_op)
        * (x_s * torch.sin(phi_op) - y_s * torch.cos(phi_op))
    )

    return grad


I_n = torch.rand((2, 2))
A = torch.rand((2, 2, 2), requires_grad=True)


int_forward = int_cost_function(A, I_n, 2)
int_forward.backward()
print("\n")
print(f"Automatic gradients:\n{A.grad}\n")
print(f"Symbolic Grad:\n{int_model_gradients(A, I_n, 2)}\n")

# Finally, now more projections. Do I have to account for rotation of sample?
#  Or is it enough that I_n is different and corresponds to A with the current rotation? Do not think so. Can sum over different axes. That is easy.
# Or use arb_projection algorithm.
# Think it is enough to sum over the other axes. That is easy.

I_n = torch.rand((3, 2, 2))
A = torch.rand((2, 2, 2), requires_grad=True)

AD_grad = torch.zeros((2, 2, 2))
for n in range(I_n.shape[0]):
    int_forward = int_cost_function(A, I_n[n], n)
    int_forward.backward()
    AD_grad += A.grad
    A.grad.zero_()
print("\n")
SYM_grad = torch.zeros((2, 2, 2))
for n in range(I_n.shape[0]):
    SYM_grad += int_model_gradients(A, I_n[n], n)

print(f"Automatic gradients:\n{AD_grad}\n")
print(f"Symbolic Grad:\n{SYM_grad}\n")


"""
Next, want to easily validate the gradients with respect to A and B, without a lengthy implementation of all the rotations. Think this is possible.
Just as in the int-model, Each I_n may represent a different projection. The corresponding Theta values are not important?
"""

simple = False
pos_phi = 0
r1, r2 = (
    2,
    3,
)  # Only applicable if simple = False and so long we only iterative over A.shape[-1]

I_n = torch.rand((3, 2, 2))
A = torch.rand((2, 2, 2), requires_grad=True)
B = torch.rand((2, 2, 2), requires_grad=True)
Theta = torch.rand((2, 2, 2), requires_grad=True)
Phi = torch.rand((2, 2, 2), requires_grad=True)

AD_grad_A = torch.zeros((2, 2, 2))
AD_grad_B = torch.zeros((2, 2, 2))
AD_grad_theta = torch.zeros((2, 2, 2))
AD_grad_phi = torch.zeros((2, 2, 2))
for n in range(r1, r2):  # I_n.shape[0]):
    forward = new_cost_function(
        A, B, Theta, Phi, I_n[n], n, simple=simple, pos_phi=pos_phi
    )
    forward.backward()
    AD_grad_A += A.grad
    AD_grad_B += B.grad
    AD_grad_theta += Theta.grad
    AD_grad_phi += Phi.grad
    A.grad.zero_()
    B.grad.zero_()
    Theta.grad.zero_()
    Phi.grad.zero_()
print("\n")

SYM_grad_A = torch.zeros((2, 2, 2))
SYM_grad_B = torch.zeros((2, 2, 2))
SYM_grad_theta = torch.zeros((2, 2, 2))
SYM_grad_phi = torch.zeros((2, 2, 2))
for n in range(r1, r2):  # I_n.shape[0]):
    SYM_grad_A += new_model_grad_A(
        A, B, Theta, Phi, I_n[n], n, simple=simple, pos_phi=pos_phi
    )
    SYM_grad_B += new_model_grad_B(
        A, B, Theta, Phi, I_n[n], n, simple=simple, pos_phi=pos_phi
    )
    SYM_grad_theta += new_model_grad_theta(
        A, B, Theta, Phi, I_n[n], n, simple=simple, pos_phi=pos_phi
    )
    SYM_grad_phi += new_model_grad_phi(
        A, B, Theta, Phi, I_n[n], n, simple=simple, pos_phi=pos_phi
    )

print("\n")
print(f"Automatic gradients A:\n{AD_grad_A}\n")
print(f"Symbolic Grad A:\n{SYM_grad_A}\n")
print(f"Automatic gradients B:\n{AD_grad_B}\n")
print(f"Symbolic Grad B:\n{SYM_grad_B}\n")
print(f"Automatic gradients Theta:\n{AD_grad_theta}\n")
print(f"Symbolic Grad Theta:\n{SYM_grad_theta}\n")
print(f"Automatic gradients Phi:\n{AD_grad_phi}\n")
print(f"Symbolic Grad Phi:\n{SYM_grad_phi}\n")


# No importing of matlab packages necessary. Located in matlab environment

# eng1 = matlab.engine.start_matlab()

# max1 = eng1.SASTT_calculate_measurement_time  # Some function

# eng1.quit()

# print(max1)  # Some output
