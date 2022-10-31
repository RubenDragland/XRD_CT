import torch
import numpy as np
import logging
import scipy.io

logger = logging.getLogger()
logger.setLevel(logging.INFO)


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
    (
        theta_struct_it,
        phi_struct_it,
        a_temp_it,
        Ylm_coef,
        ny,
        nx,
        nz,
        numOfsegments,
        numOfpixels,
        numOfCoeffs,
        numOfvoxels,
    ) = format_matlab_input(
        theta_struct_it,
        phi_struct_it,
        a_temp_it,
        Ylm_coef,
        ny,
        nx,
        nz,
        numOfsegments,
        numOfpixels,
        numOfCoeffs,
        numOfvoxels,
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

    logging.debug("error_norm: {}".format(error_norm))
    logging.debug("AD_grad_coeff shape: {}".format(AD_grad_coeff.shape))

    return format_matlab_output(error_norm, AD_grad_coeff, AD_grad_theta, AD_grad_phi)
    # Probably matlab compatible.


def format_matlab_input(
    theta_struct_it,
    phi_struct_it,
    a_temp_it,
    Ylm_coef,
    ny,
    nx,
    nz,
    numOfsegments,
    numOfpixels,
    numOfCoeffs,
    numOfvoxels,
):
    """
    Deals with matlab engine things.
    After this, assume everything is in torch.
    """
    logging.debug("Before a_temp_it shape: {}".format(a_temp_it.size))
    # RSD: Unsure on effect of reshape. Might be able to remove. Yes can be removed.
    theta_struct_it = torch.tensor(theta_struct_it.reshape(theta_struct_it.size))
    phi_struct_it = torch.tensor(phi_struct_it.reshape(phi_struct_it.size))
    a_temp_it = torch.tensor(a_temp_it.reshape(a_temp_it.size))
    Ylm_coef = torch.tensor(Ylm_coef.reshape(Ylm_coef.size))

    logging.debug("After a_temp_it shape: {}".format(a_temp_it.shape))

    # Add more transforms as needed.
    # Add possible all conversion of numbers

    ny = int(np.squeeze(ny))
    nx = int(np.squeeze(nx))
    nz = int(np.squeeze(nz))
    numOfsegments = int(np.squeeze(numOfsegments))
    numOfpixels = int(np.squeeze(numOfpixels))
    numOfCoeffs = int(np.squeeze(numOfCoeffs))
    numOfvoxels = int(np.squeeze(numOfvoxels))

    return (
        theta_struct_it,
        phi_struct_it,
        a_temp_it,
        Ylm_coef,
        ny,
        nx,
        nz,
        numOfsegments,
        numOfpixels,
        numOfCoeffs,
        numOfvoxels,
    )


def format_matlab_output(error_norm, AD_grad_coeff, AD_grad_theta, AD_grad_phi):
    """
    Prepares results for matlab engine.
    """

    error_norm = np.float64(error_norm)
    # RSD: Hopefully results in a 3D numpy array with proper size.
    if AD_grad_coeff is not None:
        AD_grad_coeff = AD_grad_coeff.cpu().detach().numpy()
    else:
        AD_grad_coeff = np.array([])
    if AD_grad_theta is not None:
        AD_grad_theta = AD_grad_theta.cpu().detach().numpy()
        AD_grad_phi = AD_grad_phi.cpu().detach().numpy()
    else:
        AD_grad_theta = np.array([])
        AD_grad_phi = np.array([])

    return (
        error_norm,
        AD_grad_coeff,
        AD_grad_theta,
        AD_grad_phi,
    )


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

    # RSD: Might have to do to.device for all tensors. data, etc. Current assumption is that it is most beneficial to only upload the tensors that are differentiated

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
    # RSD: Want data as 3D tensor. Therefore no reshape.
    data = torch.from_numpy(np.array(current_projection["data"]))

    # Rot_exp_now = np.array(
    #     current_projection["Rot_exp"].reshape(current_projection["Rot_exp"].size)
    # )  # RSD: Cannot find shape...
    Rot_exp_now = current_projection["Rot_exp"]
    unit_q_object = torch.tensor(Rot_exp_now @ unit_q_beamline)

    sin_theta_struct = torch.reshape(torch.sin(theta_struct), (1, 1, numOfvoxels))
    cos_theta_struct = torch.reshape(torch.cos(theta_struct), (1, 1, numOfvoxels))
    sin_phi_struct = torch.reshape(torch.sin(phi_struct), (1, 1, numOfvoxels))
    cos_phi_struct = torch.reshape(torch.cos(phi_struct), (1, 1, numOfvoxels))

    zeros_struct = torch.zeros((1, 1, numOfvoxels))
    ones_struct = torch.ones((1, 1, numOfvoxels))

    # Rot_str = torch.tensor(
    #     [
    #         [
    #             cos_theta_struct * cos_phi_struct,
    #             cos_theta_struct * sin_phi_struct,
    #             -sin_theta_struct,
    #         ],
    #         [-sin_phi_struct, cos_phi_struct, zeros_struct],
    #         [
    #             sin_theta_struct * cos_phi_struct,
    #             sin_theta_struct * sin_phi_struct,
    #             cos_theta_struct,
    #         ],
    #     ]
    # )
    Rot_str = torch.stack(
        [
            cos_theta_struct * cos_phi_struct,
            cos_theta_struct * sin_phi_struct,
            -sin_theta_struct,
            -sin_phi_struct,
            cos_phi_struct,
            zeros_struct,
            sin_theta_struct * cos_phi_struct,
            sin_theta_struct * sin_phi_struct,
            cos_theta_struct,
        ]
    ).reshape((3, 3, numOfvoxels))

    logging.debug("Rot_str shape: {}".format(Rot_str.shape))
    logging.debug("unit_q_object shape: {}".format(unit_q_object.shape))
    # Rot_str = torch.reshape(Rot_str, (3, 3, numOfvoxels))  # Should it be like this?

    q_pp = page_multiply(Rot_str, unit_q_object)
    # RSD: Old did not work: torch.matmul( #Rot_str, unit_q_object #)
    # RSD: Retrieve last index of q_pp
    cos_theta_sh_cut = q_pp[-1, :, :]
    logging.debug("cos_theta_sh_cut shape: {}".format(cos_theta_sh_cut.shape))

    block_cos_theta_powers = repmat_cumprod_SH(
        ones_struct, cos_theta_sh_cut, numOfCoeffs
    )

    logging.debug("Ylm_coef: {}".format(Ylm_coef))
    logging.debug(
        "block_cos_theta_powers shape: {}".format(block_cos_theta_powers.shape)
    )
    Ylm = page_multiply(Ylm_coef, block_cos_theta_powers)
    logging.debug("Ylm shape: {}".format(Ylm.shape))

    logging.debug("a_temp shape: {}".format(a_temp.shape))
    sumlm_alm_Ylm = page_multiply(a_temp, Ylm)
    logging.debug("sumlm_alm_Ylm shape: {}".format(sumlm_alm_Ylm.shape))

    # RSD: Python vs Matlab indexing
    data_synt_vol = torch.permute(torch.abs(sumlm_alm_Ylm**2), (2, 1, 0))
    data_synt_vol = torch.reshape(data_synt_vol, (ny, nx, nz, numOfsegments))

    logging.debug("data_synt_vol shape: {}".format(data_synt_vol.shape))
    logging.debug("Data type, shape: {}, {}".format(data.dtype, data.shape))

    xout = (
        np.arange(1, data.size(1) + 1)
        - np.ceil(data.size(1) / 2)
        + np.squeeze(current_projection["dx"])  # RSD: or [0]
    )
    # RSD: NB MATLAB vs Python indices
    yout = (
        np.arange(1, data.size(0) + 1)
        - np.ceil(data.size(0) / 2)
        + np.squeeze(current_projection["dy"])  # [0]
    )

    logging.debug("xout shape: {}".format(xout.shape))
    logging.debug("yout shape: {}".format(yout.shape))
    logging.debug(f"\nxout: {xout}\nyout: {yout}")

    proj_out_all = arb_projection(data_synt_vol, X, Y, Z, Rot_exp_now, p, xout, yout)

    logging.debug("data shape: {}".format(data.shape))

    permuted_data = torch.permute(data, (2, 0, 1))
    aux_diff_poisson = (
        torch.sqrt(proj_out_all) - torch.sqrt(permuted_data)
    ) * torch.from_numpy(np.array(current_projection["window_mask"]))

    logging.debug(f"aux_diff_poisson shape: {aux_diff_poisson.shape}")
    # RSD: Permute aux_diff_poisson
    aux_diff_poisson = torch.permute(
        aux_diff_poisson, (1, 2, 0)
    )  # RSD: OK? Bit unnecessary given that next is a sum over dimensions.
    error_norm = (
        2 * torch.sum(aux_diff_poisson**2, dim=(2, 1, 0)) / numOfpixels
    )  # Check dims. Dim necessary here?

    error_norm.backward()  # Backpropate gradient

    if find_coefficients:
        AD_grad_coeff = a_temp.grad
        # AD_grad_coeff = torch.autograd.grad(
        #     error_norm, a_temp, retain_graph=False
        # )  # Consider to batchify this code despite loosing parloop.
    else:
        AD_grad_coeff = None
    if find_orientation:
        AD_grad_theta = theta_struct.grad
        AD_grad_phi = phi_struct.grad
        # AD_grad_theta = torch.autograd.grad(
        #     error_norm, theta_struct, retain_graph=False
        # )
        # AD_grad_phi = torch.autograd.grad(error_norm, phi_struct, retain_graph=False)
    else:
        AD_grad_theta = None
        AD_grad_phi = None  # Might need to return zeros instead of None

    logging.debug(f"Ad_grad_coeff: {AD_grad_coeff}")
    return error_norm, AD_grad_coeff, AD_grad_theta, AD_grad_phi
    # RSD: What about shape?


def page_multiply(A: torch.tensor, B: torch.tensor):
    # Use Einstein sum in different cases
    # Needs improvements for batch, orientation, and more coefficients

    # assert len(A.size()) == 3, "A is not a 3D tensor"
    # assert len(B.size()) <= len(A.size()), "B is not a 3D or 2D tensor"

    logging.debug("A shape: {}".format(A.shape))
    logging.debug("B shape: {}".format(B.shape))

    if A.shape[-1] == B.shape[-1] and len(B.shape) == 3:

        if len(A.shape) == 1:
            A = A.unsqueeze(0).unsqueeze(0)  # RSD: Add two dimensions
        C = torch.einsum("ijm,jkm->ikm", A, B)

    elif len(A.shape) > len(B.shape):

        C = torch.einsum("ijk,jm->imk", A, B)

    elif len(A.shape) < len(B.shape) and len(B.shape) == 3:

        if len(A.shape) == 1:
            A = A.unsqueeze(0)  # RSD: Add batch dimension

        logging.debug("A dims: {}".format(len(A.shape)))
        logging.debug("B dims: {}".format(len(B.shape)))

        C = torch.einsum("ij,jkm->ikm", A, B)  # RSD: Works for a_temp as well?
    else:
        raise ValueError(
            "A and B have incompatible shapes. Please implement desired einsum"
        )
    return C


def repmat_cumprod_SH(ones_struct, cos_theta_sh_cut, numOfCoeffs):
    # RSD: This is a bit of a hack. Should be possible to do it with einsum. AI says...
    copy_matrix = torch.ones(
        (numOfCoeffs, cos_theta_sh_cut.shape[-2], cos_theta_sh_cut.shape[-1])
    )
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

    proj_out_all = torch.zeros(
        len(yout), len(xout), tomo_obj_all.shape[-1], requires_grad=False
    )
    # RSD: Check if grad must be True or False.

    logging.debug("proj_out_all shape: {}".format(proj_out_all.shape))
    logging.debug("tomo_obj_all shape: {}".format(tomo_obj_all.shape))
    nRows, nCols, nPages = proj_out_all.size()
    page_out = nRows * nCols
    page_in = Ax.shape[0] * np.product(Ax.shape[1:])
    nElements = np.product(Ax.shape)

    logging.debug(f"nRows: {nRows} nCols: {nCols} nPages: {nPages}")

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
        nCols,
        nRows,
        nPages,
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
        filter_2D = torch.unsqueeze(torch.unsqueeze(filter_2D, 0), 0)
        filter_2D = filter_2D.expand((nPages, nPages, -1, -1))

        # proj_out_all = torch.permute(
        #     torch.nn.functional.conv2d(
        #         torch.permute(proj_out_all, (2, 0, 1)), filter_2D, padding="same"
        #     ),
        #     (1, 2, 0),
        # )  # RSD: Check shape.
        proj_out_all = torch.nn.functional.conv2d(
            torch.permute(proj_out_all, (2, 0, 1)), filter_2D, padding="same"
        )
        logging.debug("Proj Shape after conv: ", proj_out_all.shape)

        # for i in range(proj_out_all.size(2)):
        #     proj_out_all[:, :, i] = torch.fft.ifft2(
        #         torch.fft.fft2(proj_out_all[:, :, i]) * torch.fft.fft2(filter_2D)
        #     )  # Convolutional theorem. Check if works. Size issue.

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

        # RSD: Matlab implementation. Will not work
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

    logging.debug("Ax shape: {}".format(Ax.shape))
    # RSD: Assume Ax is numpy array.
    indices = np.argwhere((Ax > 0) & (Ax < nCols) & (Ay > 0) & (Ay < nRows))
    arg_0, arg_1, arg_2 = indices[:, 0], indices[:, 1], indices[:, 2]

    ind = Ay[arg_0, arg_1, arg_2] + nRows * Ax[arg_0, arg_1, arg_2]

    temp1 = torch.from_numpy(Tx[arg_0, arg_1, arg_2] * Ty[arg_0, arg_1, arg_2])
    temp2 = torch.from_numpy(Tx[arg_0, arg_1, arg_2] * (1 - Ty[arg_0, arg_1, arg_2]))
    temp3 = torch.from_numpy((1 - Tx[arg_0, arg_1, arg_2]) * Ty[arg_0, arg_1, arg_2])
    temp4 = torch.from_numpy(
        (1 - Tx[arg_0, arg_1, arg_2]) * (1 - Ty[arg_0, arg_1, arg_2])
    )

    logging.debug("indices shape: {}".format(indices.shape))
    logging.debug("tomo_obj_all shape: {}".format(tomo_obj_all.shape))
    logging.debug(f"\nAy: {Ay[arg_0, arg_1, arg_2]}\nAx: {Ax[arg_0, arg_1, arg_2]}")
    logging.debug(f"temp1 shape: {temp1.shape}")
    logging.debug(
        f"tomo_obj_all indexed shape: {tomo_obj_all[arg_0, arg_1, arg_2, :].shape}"
    )

    # RSD: Believe the indexing has to be this simple
    proj_out_all[Ay[arg_0, arg_1, arg_2], Ax[arg_0, arg_1, arg_2]] += (
        tomo_obj_all[arg_0, arg_1, arg_2].T * temp1
    ).T

    proj_out_all[Ay[arg_0, arg_1, arg_2] - 1, Ax[arg_0, arg_1, arg_2]] += (
        tomo_obj_all[arg_0, arg_1, arg_2, :].T * temp2
    ).T
    proj_out_all[Ay[arg_0, arg_1, arg_2], Ax[arg_0, arg_1, arg_2] - 1] += (
        tomo_obj_all[arg_0, arg_1, arg_2, :].T * temp3
    ).T
    proj_out_all[Ay[arg_0, arg_1, arg_2] - 1, Ax[arg_0, arg_1, arg_2] - 1] += (
        tomo_obj_all[arg_0, arg_1, arg_2, :].T * temp4
    ).T
    # RSD: Here one could consider using sparse representation of the projection matrix. However, assumed to be dense enough.

    # RSD: Rage about linear indexing in MATLAB: https://www.mathworks.com/matlabcentral/answers/1015-why-does-matlab-use-linear-indexing
    # RSD: Can index as a single COLUMN vector. Hence linearly subtract one index, is to subtract one row.

    # RSD: Believe there is no need for this loop with array operations and 3D tensors. MATLAB iterates over all elements.
    """
    for j in range(nPages):
        shift_in = indices + j * page_in
        shift_out = ind + j * page_out

        proj_out_all[shift_out] += tomo_obj_all[shift_in] * temp1
        proj_out_all[shift_out - 1] += tomo_obj_all[shift_in] * temp2
        proj_out_all[shift_out - nRows] += tomo_obj_all[shift_in] * temp3
        proj_out_all[shift_out - nRows - 1] += tomo_obj_all[shift_in] * temp4
    """

    return proj_out_all


if __name__ == "__main__":
    # Test code
    workspace = scipy.io.loadmat(
        r"C:\Users\Bruker\OneDrive\Dokumenter\NTNU\XRD_CT/Data Sets/Debug Data/workspace_calc_grad.mat"
    )
    theta_struct_it = workspace["theta_struct"]
    phi_struct_it = workspace["phi_struct"]
    a_temp_it = workspace["a_temp"]
    ny = workspace["ny"]
    nx = workspace["nx"]
    nz = workspace["nz"]
    numOfsegments = workspace["numOfsegments"]
    current_projection = workspace["projection"][0, 0]  # [2]
    p = workspace["p"]
    X = workspace["X"]
    Y = workspace["Y"]
    Z = workspace["Z"]
    numOfpixels = workspace["numOfpixels"]
    unit_q_beamline = workspace["unit_q_beamline"]
    Ylm_coef = workspace["Ylm_coef"]
    find_coefficients = workspace["find_coefficients"]
    find_orientation = workspace["find_orientation"]
    numOfCoeffs = workspace["numOfCoeffs"]
    numOfvoxels = workspace["numOfvoxels"]

    main(
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
