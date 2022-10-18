import torch
import numpy as np


def print_test(string):
    print(string)
    return (string, 1, 2)


def main(
    theta_struct_it,
    phi_struct_it,
    a_temp_it,
    ny,
    nx,
    nz,
    numOfsegments,
    current_projection,
    p,
    X,
    Y,
    Z,
    numOfpixels,
    unit_q_beamline,
    Ylm_coef,
    find_coefficients,
    find_orientation,
    numOfCoeffs,
    numOfvoxels,
):
    """
    Main function.
    """
    theta_struct_it, phi_struct_it, a_temp_it = format_matlab_input(
        theta_struct_it, phi_struct_it, a_temp_it
    )

    error_norm, AD_grad_coeff, AD_grad_theta, AD_grad_phi = SAXSTT_forward_backward(
        theta_struct_it,
        phi_struct_it,
        a_temp_it,
        ny,
        nx,
        nz,
        numOfsegments,
        current_projection,
        p,
        X,
        Y,
        Z,
        numOfpixels,
        unit_q_beamline,
        Ylm_coef,
        find_coefficients,
        find_orientation,
        numOfCoeffs,
        numOfvoxels,
    )

    # Format back function here.
    error_norm = np.float64(error_norm)
    AD_grad_coeff = AD_grad_coeff.cpu().detach().numpy()
    AD_grad_theta = AD_grad_theta.cpu().detach().numpy()
    AD_grad_phi = AD_grad_phi.cpu().detach().numpy()

    return (
        error_norm,
        AD_grad_coeff,
        AD_grad_theta,
        AD_grad_phi,
    )
    # Probably matlab compatible.


def format_matlab_input(
    theta_struct_it, phi_struct_it, a_temp_it
):  # , ny, nx, nz, numOfsegments, current_projection, p, X, Y, Z, numOfpixels, unit_q_beamline, Ylm_coef, find_coefficients, find_orientation, numOfCoeffs, numOfvoxels):
    """
    Deals with matlab engine things.
    After this, assume everything is in torch.
    """

    theta_struct_it = torch.tensor(theta_struct_it.reshape(theta_struct_it.size))
    phi_struct_it = torch.tensor(phi_struct_it.reshape(phi_struct_it.size))
    a_temp_it = torch.tensor(a_temp_it.reshape(a_temp_it.size))

    # Add more transforms as needed.

    return theta_struct_it, phi_struct_it, a_temp_it


def SAXSTT_forward_backward(
    theta_struct_it,
    phi_struct_it,
    a_temp_it,
    ny,
    nx,
    nz,
    numOfsegments,
    current_projection,
    p,
    X,
    Y,
    Z,
    numOfpixels,
    unit_q_beamline,
    Ylm_coef,
    find_coefficients,
    find_orientation,
    numOfCoeffs,
    numOfvoxels,
):
    """
    Python translation of forward pass in reconstruction algorithm.
    Calculates error in reciprocal intensity given the estimated SH functions.
    The goal is to write one function in python regardless of method
    """
    if torch.cuda.is_available():
        device = torch.device("cuda")
    else:
        device = torch.device("cpu")

    if find_coefficients:
        a_temp_it.requires_grad_(True).to(device)
    if find_orientation:
        theta_struct_it.requires_grad_(True).to(device)
        phi_struct_it.requires_grad_(True).to(device)

    error_norm, AD_grad_coeff, AD_grad_theta, AD_grad_phi = SAXS_AD_cost_function(
        theta_struct_it,
        phi_struct_it,
        a_temp_it,
        ny,
        nx,
        nz,
        numOfsegments,
        current_projection,
        p,
        X,
        Y,
        Z,
        numOfpixels,
        unit_q_beamline,
        Ylm_coef,
        find_coefficients,
        find_orientation,
        numOfCoeffs,
        numOfvoxels,
    )

    return error_norm, AD_grad_coeff, AD_grad_theta, AD_grad_phi


def SAXS_AD_cost_function(
    theta_struct,
    phi_struct,
    a_temp,
    ny,
    nx,
    nz,
    numOfsegments,
    current_projection,
    p,
    X,
    Y,
    Z,
    numOfpixels,
    unit_q_beamline,
    Ylm_coef,
    find_coefficients,
    find_orientation,
    numOfCoeffs,
    numOfvoxels,
):

    data = current_projection["data"]
    Rot_exp_now = current_projection["Rot_exp_now"]
    unit_q_object = Rot_exp_now @ unit_q_beamline

    sin_theta_struct = torch.reshape(torch.sin(theta_struct), (1, 1, numOfvoxels))
    cos_theta_struct = torch.reshape(torch.cos(theta_struct), (1, 1, numOfvoxels))
    sin_phi_struct = torch.reshape(torch.sin(phi_struct), (1, 1, numOfvoxels))
    cos_phi_struct = torch.reshape(torch.cos(phi_struct), (1, 1, numOfvoxels))

    zeros_struct = torch.zeros((1, 1, numOfvoxels))
    ones_struct = torch.ones((1, 1, numOfvoxels))

    Rot_str = torch.tensor(
        [
            [
                cos_theta_struct * cos_phi_struct,
                cos_theta_struct * sin_phi_struct,
                -sin_theta_struct,
            ],
            [-sin_phi_struct, cos_phi_struct, zeros_struct],
            [
                sin_theta_struct * cos_phi_struct,
                sin_theta_struct * sin_phi_struct,
                cos_theta_struct,
            ],
        ]
    )
    Rot_str = torch.reshape(Rot_str, (3, 3, numOfvoxels))  # Should it be like this?

    q_pp = torch.matmul(Rot_str, unit_q_object)
    cos_theta_sh_cut = q_pp[3, :, :]

    block_cos_theta_powers = repmat_cumprod_SH(
        ones_struct, cos_theta_sh_cut, numOfCoeffs
    )

    Ylm = page_multiply(Ylm_coef, block_cos_theta_powers)
    sumlm_alm_Ylm = page_multiply(a_temp, Ylm)
    data_synt_vol = torch.permute(torch.abs(sumlm_alm_Ylm**2), (3, 2, 1))
    data_synt_vol = torch.reshape(data_synt_vol, (ny, nx, nz, numOfsegments))

    xout = (
        np.arange(1, data.size(1))
        - torch.ceil(data.size(1) / 2)
        + current_projection["dx"]
    )  # RSD: NB indices. Are they correct?
    yout = (
        np.arange(1, data.size(0))
        - torch.ceil(data.size(0) / 2)
        + current_projection["dy"]
    )

    proj_out_all = arb_projection(data_synt_vol, X, Y, Z, Rot_exp_now, p, xout, yout)

    aux_diff_poisson = (
        torch.sqrt(proj_out_all) - torch.sqrt(data)
    ) * current_projection["window_mask"]

    error_norm = (
        2 * torch.sum(aux_diff_poisson**2, dim=(2, 1, 0)) / numOfpixels
    )  # Check dims

    if find_coefficients:
        AD_grad_coeff = torch.autograd.grad(
            error_norm, a_temp, retain_graph=False
        )  # Consider to batchify this code despite loosing parloop.
    else:
        AD_grad_coeff = None
    if find_orientation:
        AD_grad_theta = torch.autograd.grad(
            error_norm, theta_struct, retain_graph=False
        )
        AD_grad_phi = torch.autograd.grad(error_norm, phi_struct, retain_graph=False)
    else:
        AD_grad_theta = None
        AD_grad_phi = None  # Might need to return zeros instead of None

    return error_norm, AD_grad_coeff, AD_grad_theta, AD_grad_phi
    # RSD: What about shape?


def page_multiply(A: torch.tensor, B: torch.tensor):
    # Check if numba compatible. More variables need to be converted to tensor.

    assert A.size[-1] == B.size[-1]

    for i in range(A.size[-1]):

        A[:, :, i] = torch.matmul(
            A[:, :, i], B[:, :, i]
        )  # Check if works compared to creating new Tensor.

    return A


def repmat_cumprod_SH(ones_struct, cos_theta_sh_cut, numOfCoeffs):

    copy_matrix = torch.ones((numOfCoeffs, ones_struct.size(1), ones_struct.size(2)))
    # RSD: Difference from original: Create entire matrix, iterate from 1.

    for i in range(1, numOfCoeffs):
        # RSD Edit: Start from 1, not 0. End at numOfCoeffs, not numOfCoeffs-1.
        copy_matrix[i, :, :] = copy_matrix[i, :, :] * cos_theta_sh_cut ** (2 * i)

    return copy_matrix


def arb_projection(tomo_obj_all, X, Y, Z, R, p, xout, yout):  # Assume numpy
    """
    Estimate projection based on SH. Ignore most options for now.
    """

    # Rot matrix is 3x3. X, Y, Z are meshgrid holding indices. NB! 0 indexing python vs 1 indexing matlab.
    Xp = R[0, 0] * X + R[0, 1] * Y + R[0, 2] * Z
    Yp = R[1, 0] * X + R[1, 1] * Y + R[1, 2] * Z

    # Volume upsampling ignored

    min_xout = np.min(xout)
    min_yout = np.min(yout)

    # Do bilinear interpolation directly. No need for nearest interpolation.
    Ax = np.floor(
        Xp - min_xout + 1
    )  # Drop +1 due to 0-indexing. AI suggested, but believe I should copy cpp-code
    Ay = np.floor(Yp - min_yout + 1)
    Tx = (
        Xp - min_xout + 1
    ) - Ax  # Variable from 0 to 1 from x distance of pixel Ax to Ax+1 where the voxel hits
    Ty = (Yp - min_yout + 1) - Ay

    proj_out_all = torch.zeros(len(yout), len(xout), tomo_obj_all, requires_grad=True)

    nRows, nCols, nPages = tomo_obj_all.size()
    page_out = nRows * nCols
    page_in = Ax.shape[0] * np.product(Ax.shape[1:])
    nElements = np.product(Ax.shape)

    proj_out_all = array_interpolate(
        proj_out_all,
        tomo_obj_all,
        Ax,
        Ay,
        Tx,
        Ty,
        page_in,
        page_out,
        nElements,
        nPages,
        nRows,
        nCols,
    )

    if p["filter_2D"] == 0:
        pass
    elif p["filter_2D"] == 1:
        filter_2D = torch.tensor([0.25, 1, 0.25])  # Hard code the entire matrix.
    elif p["filter_2D"] == 2:
        filter_2D = torch.tensor([0.5, 1, 0.5])
    elif p["filter_2D"] == 3:
        filter_2D = torch.tensor([1 / 3, 2 / 3, 1, 2 / 3, 1 / 3])

    if p["filter_2D"] > 0:

        filter_2D = torch.outer(filter_2D, filter_2D)

        for i in range(proj_out_all.size(2)):
            proj_out_all[:, :, i] = torch.ifft2(
                torch.fft2(proj_out_all[:, :, i]) * torch.fft2(filter_2D)
            )  # Convolutional theorem. Check if works.

    return proj_out_all


# Consider interpolation in own function.


def loop_interpolate(
    proj_out_all,
    tomo_obj_all,
    Ax,
    Ay,
    Tx,
    Ty,
    page_in,
    page_out,
    nElements,
    nCols,
    nRows,
    nPages,
):

    for i in range(
        nElements
    ):  # Try to use numba here. Or use torch functions. Or use torch functions with numba. Or array operations.

        if (Ax[i] > 0) and (Ax[i] < nCols) and (Ay[i] > 0) and (Ay[i] < nRows):
            # Check dimensions of Ax.

            ind = Ay[i] + nRows * Ax[i]
            temp1 = Tx[i] * Ty[i]
            temp2 = Tx[i] * (1 - Ty[i])
            temp3 = (1 - Tx[i]) * Ty[i]
            temp4 = (1 - Tx[i]) * (1 - Ty[i])

            for j in range(nPages):
                shift_in = i + j * page_in
                shift_out = ind + j * page_out

                proj_out_all[shift_out] += tomo_obj_all[shift_in] * temp1
                proj_out_all[shift_out - 1] += tomo_obj_all[shift_in] * temp2
                proj_out_all[shift_out - nRows] += tomo_obj_all[shift_in] * temp3
                proj_out_all[shift_out - nRows - 1] += tomo_obj_all[shift_in] * temp4

    return proj_out_all


def array_interpolate(
    proj_out_all,
    tomo_obj_all,
    Ax,
    Ay,
    Tx,
    Ty,
    page_in,
    page_out,
    nElements,
    nCols,
    nRows,
    nPages,
):
    """Interpolation with array operations to avoid for loop."""
    indices = np.argwhere((Ax > 0) and (Ax < nCols) and (Ay > 0) and (Ay < nRows))

    ind = Ay[indices] + nRows * Ax[indices]
    temp1 = Tx[indices] * Ty[indices]
    temp2 = Tx[indices] * (1 - Ty[indices])
    temp3 = (1 - Tx[indices]) * Ty[indices]
    temp4 = (1 - Tx[indices]) * (1 - Ty[indices])

    for j in range(nPages):
        shift_in = indices + j * page_in
        shift_out = ind + j * page_out

        proj_out_all[shift_out] += tomo_obj_all[shift_in] * temp1
        proj_out_all[shift_out - 1] += tomo_obj_all[shift_in] * temp2
        proj_out_all[shift_out - nRows] += tomo_obj_all[shift_in] * temp3
        proj_out_all[shift_out - nRows - 1] += tomo_obj_all[shift_in] * temp4

    return proj_out_all
