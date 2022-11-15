import torch
import numpy as np
import logging
import scipy.io
import time

logger = logging.getLogger()
logger.setLevel(logging.INFO)
# fh = logging.FileHandler("autodiff.log")
# fh.setLevel(logging.DEBUG)
# logger.addHandler(fh)


def print_test(string):
    print(string)
    return (string, 1, 2)


# RSD: Status is that array interpolation is successful, resulting in a speed-up of 60x. However, the array interpolation follows C-order, which is why the result has to be permuted.
# RSD: Should be working now.


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
    find_grad,
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
        find_grad,
    )

    logging.debug("error_norm: {}".format(error_norm))

    return format_matlab_output(error_norm, AD_grad_coeff, AD_grad_theta, AD_grad_phi)


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

    # Add more transforms as needed.
    # Add possible all conversion of numbers

    ny = int(np.squeeze(ny))
    nx = int(np.squeeze(nx))
    nz = int(np.squeeze(nz))
    numOfsegments = int(np.squeeze(numOfsegments))
    numOfpixels = int(np.squeeze(numOfpixels))
    numOfCoeffs = int(np.squeeze(numOfCoeffs))
    numOfvoxels = int(np.squeeze(numOfvoxels))

    theta_struct_it = torch.tensor(theta_struct_it)
    phi_struct_it = torch.tensor(phi_struct_it)
    a_temp_it = torch.tensor(a_temp_it)
    Ylm_coef = reshape_fortran(torch.tensor(Ylm_coef), (numOfCoeffs, numOfCoeffs))

    logging.debug("After Ylm_coef shape: {}".format(Ylm_coef.shape))

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
    find_grad,
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

    coeffs_grad = bool(find_grad and find_coefficients)
    orientation_grad = bool(find_grad and find_orientation)

    a_temp_it.requires_grad_(coeffs_grad).to(device)
    theta_struct_it.requires_grad_(orientation_grad).to(device)
    phi_struct_it.requires_grad_(orientation_grad).to(device)
    Ylm_coef.to(device)
    torch.autograd.set_detect_anomaly(True)

    # RSD: Have all tensors run on GPU if CUDA. Otherwise, CPU.

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
    device = theta_struct.device
    # RSD: Want data as 3D tensor. Therefore no reshape.
    data = torch.from_numpy(np.array(current_projection["data"])).to(device)

    Rot_exp_now = np.array(current_projection["Rot_exp"])
    unit_q_object = torch.tensor(Rot_exp_now.T @ unit_q_beamline).to(device)

    # RSD: Consider F-order reshape. Should not be necessary since it is a 1D array. However, the original array is 3D, not 1D!
    # sin_theta_struct = torch.reshape(torch.sin(theta_struct), (1, 1, numOfvoxels))
    # cos_theta_struct = torch.reshape(torch.cos(theta_struct), (1, 1, numOfvoxels))
    # sin_phi_struct = torch.reshape(torch.sin(phi_struct), (1, 1, numOfvoxels))
    # cos_phi_struct = torch.reshape(torch.cos(phi_struct), (1, 1, numOfvoxels))

    sin_theta_struct = reshape_fortran(torch.sin(theta_struct), (1, 1, numOfvoxels))
    cos_theta_struct = reshape_fortran(torch.cos(theta_struct), (1, 1, numOfvoxels))
    sin_phi_struct = reshape_fortran(torch.sin(phi_struct), (1, 1, numOfvoxels))
    cos_phi_struct = reshape_fortran(torch.cos(phi_struct), (1, 1, numOfvoxels))

    zeros_struct = torch.zeros((1, 1, numOfvoxels)).to(device)
    ones_struct = torch.ones((1, 1, numOfvoxels)).to(device)

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

    # RSD: The only case where C-order reshaping is correct.
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

    # Rot_str = reshape_fortran(Rot_str, (3, 3, numOfvoxels))

    logging.debug("Rot_str shape: {}".format(Rot_str.shape))
    logging.debug("unit_q_object shape: {}".format(unit_q_object.shape))

    q_pp = page_multiply(Rot_str, unit_q_object)

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
    data_synt_vol = reshape_fortran(data_synt_vol, (ny, nx, nz, numOfsegments))

    logging.debug("data_synt_vol shape: {}".format(data_synt_vol.shape))
    logging.debug("Data type, shape: {}, {}".format(data.dtype, data.shape))

    xout = (
        np.arange(1, data.size(1) + 1)
        - np.ceil(data.size(1) / 2)
        + np.squeeze(current_projection["dx"])
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

    aux_diff_poisson = torch.permute(
        torch.sqrt(proj_out_all) - torch.sqrt(data), (2, 0, 1)
    ) * torch.from_numpy(np.array(current_projection["window_mask"])).to(device)

    aux_diff_poisson = torch.permute(aux_diff_poisson, (1, 2, 0))
    # RSD: OK? Bit unnecessary given that next is a sum over dimensions. Gives correct error but does it affect the gradients?
    error_norm = 2 * torch.sum(aux_diff_poisson**2, dim=(2, 1, 0)) / numOfpixels

    try:
        error_norm.backward()  # Backpropate gradient
    except:
        AD_grad_coeff = None
        AD_grad_theta = None
        AD_grad_phi = None
    else:

        if find_coefficients:
            AD_grad_coeff = a_temp.grad
        else:
            AD_grad_coeff = None
        if find_orientation:
            AD_grad_theta = theta_struct.grad
            AD_grad_phi = phi_struct.grad
        else:
            AD_grad_theta = None
            AD_grad_phi = None

    finally:
        return error_norm, AD_grad_coeff, AD_grad_theta, AD_grad_phi


def page_multiply(A: torch.tensor, B: torch.tensor):
    # Use Einstein sum in different cases
    # Needs improvements for batch, orientation, and more coefficients

    logging.debug("A shape: {}".format(A.shape))
    logging.debug("B shape: {}".format(B.shape))

    if A.shape[-1] == B.shape[-1] and len(B.shape) == 3:
        # RSD: Eg sumlm_alm_Ylm
        if len(A.shape) == 1:
            A = A.unsqueeze(0).unsqueeze(0)  # RSD: Add two dimensions
        C = torch.einsum("ijm,jkm->ikm", A, B)

    elif len(A.shape) > len(B.shape):
        # RSD: Eg q_pp

        C = torch.einsum("ijk,jm->imk", A, B)

    elif len(A.shape) < len(B.shape) and len(B.shape) == 3:
        # RSD: Eg Ylm

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


# RSD: Inplace operation refuses gradients
def repmat_cumprod_SH(ones_struct, cos_theta_sh_cut, numOfCoeffs):

    device = ones_struct.device
    copy_matrix = torch.ones(
        (numOfCoeffs, cos_theta_sh_cut.shape[-2], cos_theta_sh_cut.shape[-1]),
        dtype=torch.float64,
    ).to(device)
    # RSD: Difference from original: Create entire matrix, iterate from 1.

    for i in range(1, numOfCoeffs):
        # RSD Edit: Start from 1, not 0. End at numOfCoeffs, not numOfCoeffs-1.
        copy_matrix[i, :, :] = cos_theta_sh_cut ** (2 * i)
        # copy_matrix[i, :, :] * cos_theta_sh_cut ** (2 * i)

    return copy_matrix


def arb_projection(tomo_obj_all, X, Y, Z, R, p, xout, yout):  # Assume numpy
    """
    Estimate projection based on SH. Ignore most options for now.
    """
    device = tomo_obj_all.device
    # Rot matrix is 3x3. X, Y, Z are meshgrid holding indices. NB! 0 indexing python vs 1 indexing matlab.
    Xp = R[0, 0] * X + R[0, 1] * Y + R[0, 2] * Z
    Yp = R[1, 0] * X + R[1, 1] * Y + R[1, 2] * Z

    # Volume upsampling ignored

    min_xout = np.min(xout)
    min_yout = np.min(yout)

    # Do bilinear interpolation directly. No need for nearest interpolation.

    # Drop +1 due to 0-indexing. AI suggested, but believe I should copy cpp-code
    Ax = np.floor(Xp - min_xout + 1)
    Ay = np.floor(Yp - min_yout + 1)
    # Variable from 0 to 1 from x distance of pixel Ax to Ax+1 where the voxel hits
    Tx = (Xp - min_xout + 1) - Ax
    Ty = (Yp - min_yout + 1) - Ay

    proj_out_all = torch.zeros(
        len(yout),
        len(xout),
        tomo_obj_all.shape[-1],
        requires_grad=False,
        dtype=torch.float64,
    ).to(device)

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
        filter_2D = torch.tensor([0.25, 1, 0.25], dtype=torch.float64).to(device)
    elif p["filter_2D"] == 2:
        filter_2D = torch.tensor([0.5, 1, 0.5], dtype=torch.float64).to(device)
    elif p["filter_2D"] == 3:
        filter_2D = torch.tensor(
            [1 / 3, 2 / 3, 1, 2 / 3, 1 / 3], dtype=torch.float64
        ).to(device)

    if p["filter_2D"] > 0:

        filter_2D = torch.outer(filter_2D, filter_2D)
        filter_2D /= torch.sum(filter_2D)
        logging.debug(f"filter_2D: {filter_2D}")
        filter_2D = torch.unsqueeze(torch.unsqueeze(filter_2D, 0), 0)
        filter_2D = filter_2D.expand((nPages, -1, -1, -1))
        logging.debug(f"filter_2D shape: {filter_2D.shape}")
        logging.debug(f"full filter 2D: {filter_2D}")
        logging.debug(f"proj_out_all before conv shape: {proj_out_all.shape}")

        proj_out_all = torch.unsqueeze(proj_out_all, 0)
        proj_out_all = torch.permute(proj_out_all, (0, 3, 1, 2))
        logging.debug(f"proj_out_all ready conv: {proj_out_all[0,0,:,:]}")
        proj_out_all = torch.nn.functional.conv2d(
            proj_out_all, filter_2D, padding="same", groups=nPages
        )
        proj_out_all = torch.squeeze(torch.permute(proj_out_all, (0, 2, 3, 1)), 0)
        logging.debug(f"proj_out_all after conv shape: {proj_out_all.shape}")

    return proj_out_all


# Slow but works.
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
    # RSD: Issue here to replicate MATLAB shapes. It indexes along the first column first, not the first row. After all columns, begins the third dimension (I Think)
    # RSD: Fixed with borrowed Fortran reshape from STACKOVERFLOW. Numpy has this ability supported.
    org_shape = proj_out_all.shape
    proj_out_all = reshape_fortran(
        proj_out_all, (-1,)
    )  # torch.reshape(proj_out_all, (-1,))
    tomo_obj_all = reshape_fortran(
        tomo_obj_all, (-1,)
    )  # torch.reshape(tomo_obj_all, (-1,))
    Ax = np.reshape(Ax, (-1,), order="F")
    Ay = np.reshape(Ay, (-1,), order="F")
    Tx = np.reshape(Tx, (-1,), order="F")
    Ty = np.reshape(Ty, (-1,), order="F")

    for i in range(
        nElements
    ):  # Try to use numba here. Or use torch functions. Or use torch functions with numba. Or array operations.

        # RSD: Matlab implementation. Will not work
        if (Ax[i] > 0) & (Ax[i] < nCols) & (Ay[i] > 0) & (Ay[i] < nRows):
            # Check dimensions of Ax.
            ind = Ay[i] + nRows * Ax[i]
            temp1 = Tx[i] * Ty[i]
            temp2 = Tx[i] * (1 - Ty[i])
            temp3 = (1 - Tx[i]) * Ty[i]
            temp4 = (1 - Tx[i]) * (1 - Ty[i])

            for j in range(nPages):
                shift_in = int(i + j * page_in)
                shift_out = int(ind + j * page_out)

                proj_out_all[shift_out] = (
                    proj_out_all[shift_out] + tomo_obj_all[shift_in] * temp1
                )
                proj_out_all[shift_out - 1] = (
                    proj_out_all[shift_out - 1] + tomo_obj_all[shift_in] * temp2
                )
                proj_out_all[shift_out - nRows] = (
                    proj_out_all[shift_out - nRows] + tomo_obj_all[shift_in] * temp3
                )
                proj_out_all[shift_out - nRows - 1] = (
                    proj_out_all[shift_out - nRows - 1] + tomo_obj_all[shift_in] * temp4
                )

                # print(tomo_obj_all[shift_in] * temp4)

    proj_out_all = reshape_fortran(proj_out_all, org_shape)
    # proj_out_all.reshape(org_shape, order="F")
    logging.debug(f"proj_out_all before conv:\n{proj_out_all}")
    return proj_out_all


# RSD: Fast and functional
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

    device = proj_out_all.device
    logging.debug("Ax shape: {}".format(Ax.shape))

    indices = np.argwhere((Ax > 0) & (Ax < nCols) & (Ay > 0) & (Ay < nRows))
    arg_0, arg_1, arg_2 = indices[:, 0], indices[:, 1], indices[:, 2]

    temp1 = (
        torch.from_numpy(Tx[arg_0, arg_1, arg_2] * Ty[arg_0, arg_1, arg_2]).unsqueeze(1)
    ).to(device)
    temp2 = (
        torch.from_numpy(
            Tx[arg_0, arg_1, arg_2] * (1 - Ty[arg_0, arg_1, arg_2])
        ).unsqueeze(1)
    ).to(device)
    temp3 = (
        torch.from_numpy(
            (1 - Tx[arg_0, arg_1, arg_2]) * Ty[arg_0, arg_1, arg_2]
        ).unsqueeze(1)
    ).to(device)
    temp4 = (
        torch.from_numpy(
            (1 - Tx[arg_0, arg_1, arg_2]) * (1 - Ty[arg_0, arg_1, arg_2])
        ).unsqueeze(1)
    ).to(device)

    logging.debug("indices shape: {}".format(indices.shape))
    logging.debug("tomo_obj_all shape: {}".format(tomo_obj_all.shape))
    logging.debug(f"temp1 shape: {temp1[0,:]}")
    logging.debug(
        f"tomo_obj_all indexed shape: {tomo_obj_all[arg_1, arg_0, arg_2, :].shape}"
    )
    # RSD: The issue has been that indexing is occuring in parallel, not accumulative. Hence, a special formula is needed.
    Ax = torch.from_numpy(Ax.astype(np.int64)).to(device)
    Ay = torch.from_numpy(Ay.astype(np.int64)).to(device)

    proj_out_all = torch.index_put(
        proj_out_all,
        (Ay[arg_0, arg_1, arg_2], Ax[arg_0, arg_1, arg_2]),
        tomo_obj_all[arg_0, arg_1, arg_2, :] * temp1,
        accumulate=True,
    )
    proj_out_all = torch.index_put(
        proj_out_all,
        (Ay[arg_0, arg_1, arg_2] - 1, Ax[arg_0, arg_1, arg_2]),
        tomo_obj_all[arg_0, arg_1, arg_2, :] * temp2,
        accumulate=True,
    )
    proj_out_all = torch.index_put(
        proj_out_all,
        (Ay[arg_0, arg_1, arg_2], Ax[arg_0, arg_1, arg_2] - 1),
        tomo_obj_all[arg_0, arg_1, arg_2, :] * temp3,
        accumulate=True,
    )
    proj_out_all = torch.index_put(
        proj_out_all,
        (Ay[arg_0, arg_1, arg_2] - 1, Ax[arg_0, arg_1, arg_2] - 1),
        tomo_obj_all[arg_0, arg_1, arg_2, :] * temp4,
        accumulate=True,
    )
    # RSD: Rage about linear indexing in MATLAB: https://www.mathworks.com/matlabcentral/answers/1015-why-does-matlab-use-linear-indexing
    # RSD: Can index as a single COLUMN vector. Hence linearly subtract one index, is to subtract one row.

    return proj_out_all


def reshape_fortran(x, shape):
    if len(x.shape) > 0:
        x = x.permute(*reversed(range(len(x.shape))))
    return x.reshape(*reversed(shape)).permute(*reversed(range(len(shape))))


if __name__ == "__main__":
    # Test code
    workspace = scipy.io.loadmat(
        r"C:\Users\Bruker\OneDrive\Dokumenter\NTNU\XRD_CT\Data sets\Debug Data\orientation_python_workspace.mat"
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

    tic = time.time()

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
        find_grad=True,
    )
    tac = time.time()
    print(f"Time elapsed: {tac - tic}")

# RSD: False code:
"""
    # RSD: Believe the indexing has to be this simple. Tomo-indexing has to be wrong. tomo have to swap order of arg1 and arg0
    # proj_out_all[Ay[arg_0, arg_1, arg_2], Ax[arg_0, arg_1, arg_2], :] = (
    #     proj_out_all[Ay[arg_0, arg_1, arg_2], Ax[arg_0, arg_1, arg_2], :]
    #     + temp1 * tomo_obj_all[arg_1, arg_0, arg_2, :]
    # )

    # proj_out_all[Ay[arg_0, arg_1, arg_2] - 1, Ax[arg_0, arg_1, arg_2]] = (
    #     proj_out_all[Ay[arg_0, arg_1, arg_2] - 1, Ax[arg_0, arg_1, arg_2]]
    #     + temp2 * tomo_obj_all[arg_1, arg_0, arg_2, :]
    # )

    # print(torch.max(temp2.unsqueeze(1) * tomo_obj_all[arg_1, arg_0, arg_2, :]))
    # proj_out_all[Ay[arg_0, arg_1, arg_2], Ax[arg_0, arg_1, arg_2] - 1] = (
    #     proj_out_all[Ay[arg_0, arg_1, arg_2], Ax[arg_0, arg_1, arg_2] - 1]
    #     + temp3 * tomo_obj_all[arg_1, arg_0, arg_2, :]
    # )

    # proj_out_all[Ay[arg_0, arg_1, arg_2] - 1, Ax[arg_0, arg_1, arg_2] - 1] = (
    #     proj_out_all[Ay[arg_0, arg_1, arg_2] - 1, Ax[arg_0, arg_1, arg_2] - 1]
    #     + temp4 * tomo_obj_all[arg_1, arg_0, arg_2, :]
    # )


"""
