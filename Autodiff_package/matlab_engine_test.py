import matlab.engine
import numpy as np
import sys
import os
import torch

# Necessary?
def R_matrix(alpha):
    Rn_exp = torch.tensor(
        [
            [torch.cos(alpha), 0, torch.sin(alpha)],
            [0, 1, 0],
            [-torch.sin(alpha), 0, torch.cos(alpha)],
        ]
    )
    return Rn_exp


def R_matrix_inv(alpha):
    Rn_exp = torch.tensor(
        [
            [torch.cos(alpha), 0, -torch.sin(alpha)],
            [0, 1, 0],
            [torch.sin(alpha), 0, torch.cos(alpha)],
        ]
    )
    return Rn_exp


def exp_sin_squared(A, B, Theta):
    return A**2 * torch.exp(-B * torch.sin(Theta) ** 2)


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


def derive_SIN_THETA(A, B, theta_op, phi_op, n):
    # Derive the Theta function
    alpha, beta = lab_rotation_given_axis_n(n)

    # Rotation matrices
    COS_THETA = torch.sin(theta_op) * torch.cos(phi_op) * torch.cos(beta) + torch.cos(
        theta_op
    ) * torch.sin(beta)

    #
    return


def new_cost_function(A, B, theta, phi, I_n, n):
    # alpha, beta = lab_rotation_given_axis_n(n)
    Theta = theta
    # Temp rotation # TODO: Get from matlab. Implement properly... Does not necessarily make sense. But will be the same issue in both implementations
    I_n_hat = torch.sum(exp_sin_squared(A, B, Theta), axis=n)
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


def BP(A, B, Theta, I_n, n):
    return (
        torch.sum(exp_sin_squared(A, B, Theta), axis=n) ** 0.5 - I_n**0.5
    ) / torch.sum(exp_sin_squared(A, B, Theta), axis=n) ** 0.5


def new_model_grad_A(A, B, Theta, I_n, n):
    grad_A = torch.zeros(A.shape)
    BP_n = BP(A, B, Theta, I_n, n)
    for i in range(A.shape[-1]):
        grad_A[:, :, i] = BP_n * exp_sin_squared(A, B, Theta)[:, :, i] / A[:, :, i]
        # More efficient if index params. Note A**2 term included in exp_sin_squared function. Have to divide by A
    return torch.permute(
        4 * grad_A, (2, 0, 1) if n == 0 else (0, 2, 1) if n == 1 else (0, 1, 2)
    )


def new_model_grad_B(A, B, Theta, I_n, n):
    grad_B = torch.zeros(A.shape)
    BP_n = BP(A, B, Theta, I_n, n)
    for i in range(A.shape[-1]):
        grad_B[:, :, i] = (
            BP_n
            * torch.sin(Theta)[:, :, i] ** 2
            * exp_sin_squared(A, B, Theta)[:, :, i]
        )
        # More efficient if index params. Note A**2 term included in exp_sin_squared function.
    return torch.permute(
        -2 * grad_B, (2, 0, 1) if n == 0 else (0, 2, 1) if n == 1 else (0, 1, 2)
    )


# Ignore the rotation for now. Test one projection first.

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

I_n = torch.rand((3, 2, 2))
A = torch.rand((2, 2, 2), requires_grad=True)
B = torch.rand((2, 2, 2), requires_grad=True)
Theta = torch.rand((3, 2, 2, 2))

AD_grad_A = torch.zeros((2, 2, 2))
AD_grad_B = torch.zeros((2, 2, 2))
for n in range(0, 3):  # I_n.shape[0]):
    forward = new_cost_function(A, B, Theta[n], 0, I_n[n], n)
    forward.backward()
    AD_grad_A += A.grad
    AD_grad_B += B.grad
    A.grad.zero_()
    B.grad.zero_()
print("\n")

SYM_grad_A = torch.zeros((2, 2, 2))
SYM_grad_B = torch.zeros((2, 2, 2))
for n in range(0, 3):  # I_n.shape[0]):
    SYM_grad_A += new_model_grad_A(A, B, Theta[n], I_n[n], n)
    SYM_grad_B += new_model_grad_B(A, B, Theta[n], I_n[n], n)

print(f"Automatic gradients A:\n{AD_grad_A}\n")
print(f"Symbolic Grad A:\n{SYM_grad_A}\n")
print(f"Automatic gradients B:\n{AD_grad_B}\n")
print(f"Symbolic Grad B:\n{SYM_grad_B}\n")


sys.exit()


A = torch.tensor([3.0], requires_grad=True)
B = torch.tensor([0.5], requires_grad=True)
Theta = torch.tensor([0.33], requires_grad=True)


fval = exp_sin_squared(A, B, Theta)
fval.backward()
grad_A, grad_B, grad_T = A.grad, B.grad, Theta.grad

print(f"{fval}\n{grad_A}\n{grad_B}\n{grad_T}")


# No importing of matlab packages necessary. Located in matlab environment

# eng1 = matlab.engine.start_matlab()

# max1 = eng1.SASTT_calculate_measurement_time  # Some function

# eng1.quit()

# print(max1)  # Some output
