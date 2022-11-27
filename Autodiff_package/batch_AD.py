import torch
import numpy as np
import logging
import scipy.io
import time

logger = logging.getLogger()
logger.setLevel(logging.INFO)


def SH_main(
    theta_struct_it,
    phi_struct_it,
    a_temp_it,
    ny,
    nx,
    nz,
    numOfsegments,
    projections,
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
    ) = format_matlab_input_SH(
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
    if projections == 1:
        try:
            p_projection_filename = p["projection_filename"][0, 0][0]
            projections = scipy.io.loadmat(rf"{p_projection_filename}")["projection"][0]
        except:
            p_projection_filename = p["projection_filename"]
            projections = scipy.io.loadmat(rf"{p_projection_filename}")["projection"][0]
        finally:
            pass
    else:
        projections = scipy.io.loadmat(rf"{projections}")["projection"][0]

    error_norm, AD_grad_coeff, AD_grad_theta, AD_grad_phi = SH_forward_backward(
        theta_struct_it,
        phi_struct_it,
        a_temp_it,
        ny,
        nx,
        nz,
        numOfsegments,
        projections,
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

    return format_matlab_output_SH(
        error_norm, AD_grad_coeff, AD_grad_theta, AD_grad_phi
    )


def EXPSIN_main(
    theta,
    phi,
    A,
    B,
    unit_q_beamline,
    p,
    X,
    Y,
    Z,
    ny,
    nx,
    nz,
    numOfsegments,
    numOfpixels,
    numOfvoxels,
    find_grad,
):
    """
    Main function of forward and backward propagation for EXP_SIN_SQARED Cost function.
    """
    # RSD: Depending on working from python in debugging or running matlab actual reconstructions.
    try:
        p_projection_filename = p["projection_filename"][0, 0][0]
        projections = scipy.io.loadmat(rf"{p_projection_filename}")["projection"][0]
    except:
        p_projection_filename = p["projection_filename"]
        projections = scipy.io.loadmat(rf"{p_projection_filename}")["projection"][0]
    finally:
        pass

    (
        theta,
        phi,
        A,
        B,
        ny,
        nx,
        nz,
        numOfsegments,
        numOfpixels,
        numOfvoxels,
    ) = format_matlab_input_EXPSIN(
        theta,
        phi,
        A,
        B,
        ny,
        nx,
        nz,
        numOfsegments,
        numOfpixels,
        numOfvoxels,
    )

    (
        error_norm,
        AD_grad_A,
        AD_grad_B,
        AD_grad_theta,
        AD_grad_phi,
    ) = EXPSIN_forward_backward(
        theta,
        phi,
        A,
        B,
        projections,
        unit_q_beamline,
        p,
        X,
        Y,
        Z,
        ny,
        nx,
        nz,
        numOfsegments,
        numOfpixels,
        numOfvoxels,
        find_grad,
    )

    return format_matlab_output_EXPSIN(
        error_norm, AD_grad_A, AD_grad_B, AD_grad_theta, AD_grad_phi
    )


def format_matlab_input_EXPSIN(
    theta,
    phi,
    A,
    B,
    ny,
    nx,
    nz,
    numOfsegments,
    numOfpixels,
    numOfvoxels,
):

    ny = int(np.squeeze(ny))
    nx = int(np.squeeze(nx))
    nz = int(np.squeeze(nz))
    numOfsegments = int(np.squeeze(numOfsegments))
    numOfpixels = int(np.squeeze(numOfpixels))
    numOfvoxels = int(np.squeeze(numOfvoxels))

    theta = torch.tensor(theta, dtype=torch.float64)
    phi = torch.tensor(phi, dtype=torch.float64)
    A = torch.tensor(A, dtype=torch.float64)
    B = torch.tensor(B, dtype=torch.float64)

    return (
        theta,
        phi,
        A,
        B,
        ny,
        nx,
        nz,
        numOfsegments,
        numOfpixels,
        numOfvoxels,
    )


def format_matlab_input_SH(
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


# RSD: Bit lazy to copy so much code.
def format_matlab_output_EXPSIN(
    error_norm, AD_grad_A, AD_grad_B, AD_grad_theta, AD_grad_phi
):
    error_norm = np.float64(error_norm)
    if AD_grad_A is not None:
        AD_grad_A = AD_grad_A.cpu().detach().numpy()
        AD_grad_B = AD_grad_B.cpu().detach().numpy()
    else:
        AD_grad_A = np.array([])
        AD_grad_B = np.array([])
    if AD_grad_theta is not None:
        AD_grad_theta = AD_grad_theta.cpu().detach().numpy()
        AD_grad_phi = AD_grad_phi.cpu().detach().numpy()
    else:
        AD_grad_theta = np.array([])
        AD_grad_phi = np.array([])

    return (
        error_norm,
        AD_grad_A,
        AD_grad_B,
        AD_grad_theta,
        AD_grad_phi,
    )


def format_matlab_output_SH(error_norm, AD_grad_coeff, AD_grad_theta, AD_grad_phi):
    """
    Prepares results for matlab engine.
    """

    error_norm = np.float64(error_norm)
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


def reshape_fortran(x, shape):
    if len(x.shape) > 0:
        x = x.permute(*reversed(range(len(x.shape))))
    return x.reshape(*reversed(shape)).permute(*reversed(range(len(shape))))


def reshape_projections(proj, n_proj, key="data"):
    """
    Reshapes projections to be in the right format. Optimised version will not work.
    # indices = np.arange(len(proj), dtype=int)
    # projections[indices] = proj[indices] Does not seem to work
    """

    if key == "data":
        projections = np.zeros(
            (n_proj, proj[0].shape[0], proj[0].shape[1], proj[0].shape[2])
        )
    elif key == "window_mask":
        projections = np.zeros((n_proj, proj[0].shape[0], proj[0].shape[1]))
    elif key == "Rot_exp":
        projections = np.zeros((n_proj, 3, 3))
    elif key == "dx" or key == "dy":
        projections = np.zeros((n_proj, 1))
    else:
        raise ValueError("Unknown key: {}".format(key))

    for i in range(n_proj):
        projections[i] = proj[i]

    return projections


def EXPSIN_forward_backward(
    theta,
    phi,
    A,
    B,
    projections,
    unit_q_beamline,
    p,
    X,
    Y,
    Z,
    ny,
    nx,
    nz,
    numOfsegments,
    numOfpixels,
    numOfvoxels,
    find_grad,
):
    if torch.cuda.is_available():
        device = torch.device("cuda")
    else:
        device = torch.device("cpu")

    find_grad = bool(find_grad)

    theta.requires_grad_(find_grad).to(device)
    phi.requires_grad_(find_grad).to(device)
    A.requires_grad_(find_grad).to(device)
    B.requires_grad_(find_grad).to(device)
    torch.autograd.set_detect_anomaly(True)

    (
        error_norm,
        AD_grad_A,
        AD_grad_B,
        AD_grad_theta,
        AD_grad_phi,
    ) = EXPSIN_AD_cost_function(
        theta,
        phi,
        A,
        B,
        projections,
        unit_q_beamline,
        p,
        X,
        Y,
        Z,
        ny,
        nx,
        nz,
        numOfsegments,
        numOfpixels,
        numOfvoxels,
        find_grad,
    )

    return error_norm, AD_grad_A, AD_grad_B, AD_grad_theta, AD_grad_phi


def SH_forward_backward(
    theta_struct_it,
    phi_struct_it,
    a_temp_it,
    ny,
    nx,
    nz,
    numOfsegments,
    projections,
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

    error_norm, AD_grad_coeff, AD_grad_theta, AD_grad_phi = SH_AD_cost_function(
        theta_struct_it,
        phi_struct_it,
        a_temp_it,
        ny,
        nx,
        nz,
        numOfsegments,
        projections,
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


def EXPSIN_AD_cost_function(
    theta,
    phi,
    A,
    B,
    projections,
    unit_q_beamline,
    p,
    X,
    Y,
    Z,
    ny,
    nx,
    nz,
    numOfsegments,
    numOfpixels,
    numOfvoxels,
    find_grad,
):
    # RSD: This part could be in a common function...
    device = theta.device
    n_proj = len(projections)  # Hope so.
    data = (
        torch.from_numpy(
            reshape_projections(projections["data"], n_proj, key="data"),
        )
        .to(device)
        .double()
    )
    Rot_exp_now = np.array(
        reshape_projections(projections["Rot_exp"], n_proj, key="Rot_exp")
    )

    unit_q_object = torch.tensor(
        np.transpose(Rot_exp_now, (0, 2, 1)) @ unit_q_beamline
    ).to(device)

    sin_theta = reshape_fortran(torch.sin(theta), (1, 1, numOfvoxels))
    cos_theta = reshape_fortran(torch.cos(theta), (1, 1, numOfvoxels))
    sin_phi = reshape_fortran(torch.sin(phi), (1, 1, numOfvoxels))
    cos_phi = reshape_fortran(torch.cos(phi), (1, 1, numOfvoxels))

    zeros_struct = torch.zeros((1, 1, numOfvoxels)).to(device)

    Rot_str = torch.stack(
        [
            cos_theta * cos_phi,
            cos_theta * sin_phi,
            -sin_theta,
            -sin_phi,
            cos_phi,
            zeros_struct,
            sin_theta * cos_phi,
            sin_theta * sin_phi,
            cos_theta,
        ]
    ).reshape((3, 3, numOfvoxels))

    q_pp = batch_multiply(Rot_str, unit_q_object)
    COS_THETA = q_pp[:, -1, :, :]
    COS_THETA = torch.permute(COS_THETA, (0, 2, 1))
    A_sq = torch.squeeze(A).unsqueeze(0).unsqueeze(-1)
    B_sq = torch.squeeze(B).unsqueeze(0).unsqueeze(-1)
    # RSD: Have to create a new tensor to keep original variables leaf.

    # RSD: Here the new things begin:

    # RSD: Normalisation put on hold for now. Would have to be implemented manually...

    norm = torch.sqrt(torch_simpsons(norm_integrand, B_sq) ** 2)  # RSD: Check shape

    SIN_THETA_SQUARED = 1 - COS_THETA**2

    I_hat_tomogram = A_sq**2 * torch.exp(-B_sq * SIN_THETA_SQUARED).to(device) / norm
    I_hat_tomogram = reshape_fortran(
        I_hat_tomogram, (n_proj, ny, nx, nz, numOfsegments)
    )

    dx = reshape_projections(projections["dx"], n_proj, "dx")
    dy = reshape_projections(projections["dy"], n_proj, "dy")

    xout = (
        np.expand_dims(
            np.arange(1, data.size(2) + 1) - np.ceil(data.size(2) / 2), axis=0
        )
        + dx
    )
    yout = (
        np.expand_dims(
            np.arange(1, data.size(1) + 1) - np.ceil(data.size(1) / 2), axis=0
        )
        + dy
    )

    proj_out_all = arb_projection(I_hat_tomogram, X, Y, Z, Rot_exp_now, p, xout, yout)

    aux_diff_poisson = (
        torch.sqrt(proj_out_all + 1e-12) - torch.sqrt(data)
    ) * torch.from_numpy(
        reshape_projections(projections["window_mask"], n_proj, "window_mask")[
            :, :, :, np.newaxis
        ]
    ).to(
        device
    )

    aux_diff_poisson = torch.permute(aux_diff_poisson, (0, 2, 3, 1))
    error_norm = 2 * torch.sum(aux_diff_poisson**2, dim=(3, 2, 1, 0)) / numOfpixels

    if not find_grad:
        return error_norm, None, None, None, None

    try:
        error_norm.backward()  # Backpropate gradient
    except:
        AD_grad_A = None
        AD_grad_B = None
        AD_grad_theta = None
        AD_grad_phi = None
    else:
        AD_grad_A = A.grad
        AD_grad_B = B.grad
        AD_grad_theta = theta.grad
        AD_grad_phi = phi.grad

    finally:
        return error_norm, AD_grad_A, AD_grad_B, AD_grad_theta, AD_grad_phi


def SH_AD_cost_function(
    theta_struct,
    phi_struct,
    a_temp,
    ny,
    nx,
    nz,
    numOfsegments,
    projections,
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
    n_proj = len(
        projections
    )  # projections.shape[0]  # RSD: Check this. Should be 255 but might be 1
    # RSD: Want data as 3D tensor. Therefore no reshape.

    data = (
        torch.from_numpy(
            np.array(reshape_projections(projections["data"], n_proj, "data")),
        )
        .to(device)
        .double()
    )

    Rot_exp_now = np.array(
        reshape_projections(projections["Rot_exp"], n_proj, "Rot_exp")
    )

    # RSD: Shape (n_proj, 3,8)
    unit_q_object = torch.tensor(
        np.transpose(Rot_exp_now, (0, 2, 1)) @ unit_q_beamline
    ).to(device)
    # RSD: Will have to be changed Ensure n projections of that matrix multiplication. Update: Supposed to work.

    sin_theta_struct = reshape_fortran(torch.sin(theta_struct), (1, 1, numOfvoxels))
    cos_theta_struct = reshape_fortran(torch.cos(theta_struct), (1, 1, numOfvoxels))
    sin_phi_struct = reshape_fortran(torch.sin(phi_struct), (1, 1, numOfvoxels))
    cos_phi_struct = reshape_fortran(torch.cos(phi_struct), (1, 1, numOfvoxels))

    zeros_struct = torch.zeros((1, 1, numOfvoxels)).to(device)

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

    q_pp = batch_multiply(Rot_str, unit_q_object)

    # RSD: Retrieve last index of q_pp
    cos_theta_sh_cut = q_pp[:, -1, :, :]  # RSD: Shape (n_proj,...)
    logging.debug("cos_theta_sh_cut shape: {}".format(cos_theta_sh_cut.shape))

    # RSD: Shape (n_proj, numOfCoeffs, ... )
    block_cos_theta_powers = repmat_cumprod_SH(cos_theta_sh_cut, numOfCoeffs, n_proj)

    logging.debug(
        "block_cos_theta_powers shape: {}".format(block_cos_theta_powers.shape)
    )
    Ylm = batch_multiply(Ylm_coef, block_cos_theta_powers)
    logging.debug("Ylm shape: {}".format(Ylm.shape))

    logging.debug("a_temp shape: {}".format(a_temp.shape))
    sumlm_alm_Ylm = batch_multiply(a_temp, Ylm)
    logging.debug("sumlm_alm_Ylm shape: {}".format(sumlm_alm_Ylm.shape))

    # RSD: Python vs Matlab indexing. Can become challenging.
    data_synt_vol = torch.permute(torch.abs(sumlm_alm_Ylm**2), (0, 3, 2, 1))
    data_synt_vol = reshape_fortran(data_synt_vol, (n_proj, ny, nx, nz, numOfsegments))

    dx = reshape_projections(projections["dx"], n_proj, "dx")
    dy = reshape_projections(projections["dy"], n_proj, "dy")

    # Check data size. This might now be 2, 1 instead.
    xout = (
        np.expand_dims(
            np.arange(1, data.size(2) + 1) - np.ceil(data.size(2) / 2), axis=0
        )
        + dx
    )
    yout = (
        np.expand_dims(
            np.arange(1, data.size(1) + 1) - np.ceil(data.size(1) / 2), axis=0
        )
        + dy
    )

    logging.debug("xout shape: {}".format(xout.shape))
    logging.debug("yout shape: {}".format(yout.shape))
    logging.debug(f"\nxout: {xout}\nyout: {yout}")

    # RSD: Finished until here. shape n_proj, ny, nx, nz, numOfsegments
    proj_out_all = arb_projection(data_synt_vol, X, Y, Z, Rot_exp_now, p, xout, yout)

    logging.debug("data shape: {}".format(data.shape))

    aux_diff_poisson = (
        torch.sqrt(proj_out_all + 1e-12) - torch.sqrt(data)
    ) * torch.from_numpy(
        reshape_projections(projections["window_mask"], n_proj, "window_mask")[
            :, :, :, np.newaxis
        ]
    ).to(
        device
    )

    aux_diff_poisson = torch.permute(aux_diff_poisson, (0, 2, 3, 1))
    # RSD: OK? Should have nothing to say.
    error_norm = 2 * torch.sum(aux_diff_poisson**2, dim=(3, 2, 1, 0)) / numOfpixels

    try:
        error_norm.backward()  # Backpropate gradient
        pass
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


def arb_projection(tomo_obj_all, X, Y, Z, R, p, xout, yout):  # Assume numpy
    """
    Estimate projection based on SH. Ignore most options for now.
    """
    device = tomo_obj_all.device
    n_proj = tomo_obj_all.size(0)
    X, Y, Z = (
        np.expand_dims(X, axis=0),
        np.expand_dims(Y, axis=0),
        np.expand_dims(Z, axis=0),
    )
    R = R[:, :, :, np.newaxis, np.newaxis, np.newaxis]  # np.expand_dims(R, axis=3)

    # Rot matrix is 3x3. X, Y, Z are meshgrid holding indices. NB! 0 indexing python vs 1 indexing matlab.
    Xp = R[:, 0, 0] * X + R[:, 0, 1] * Y + R[:, 0, 2] * Z
    Yp = R[:, 1, 0] * X + R[:, 1, 1] * Y + R[:, 1, 2] * Z

    # RSD: One value for each projection. Believe this is along axis 1.
    min_xout = np.min(xout, axis=1)[:, np.newaxis, np.newaxis, np.newaxis]
    min_yout = np.min(yout, axis=1)[:, np.newaxis, np.newaxis, np.newaxis]

    # RSD: Think is broadcastable
    Ax = np.floor(Xp - min_xout + 1)
    Ay = np.floor(Yp - min_yout + 1)
    # Variable from 0 to 1 from x distance of pixel Ax to Ax+1 where the voxel hits
    Tx = (Xp - min_xout + 1) - Ax
    Ty = (Yp - min_yout + 1) - Ay

    proj_out_all = torch.zeros(
        R.shape[0],
        len(yout[0]),
        len(xout[0]),
        tomo_obj_all.shape[-1],
        requires_grad=False,
        dtype=torch.float64,
    ).to(device)

    logging.debug("proj_out_all shape: {}".format(proj_out_all.shape))
    logging.debug("tomo_obj_all shape: {}".format(tomo_obj_all.shape))
    n_proj, nRows, nCols, nPages = proj_out_all.size()

    logging.debug(f"nRows: {nRows} nCols: {nCols} nPages: {nPages}")

    proj_out_all = array_interpolate(
        proj_out_all,
        tomo_obj_all,
        Ax,
        Ay,
        Tx,
        Ty,
        nCols,
        nRows,
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
        filter_2D = torch.unsqueeze(torch.unsqueeze(filter_2D, 0), 0)
        filter_2D = filter_2D.expand((nPages, -1, -1, -1))
        logging.debug(f"filter_2D shape: {filter_2D.shape}")
        logging.debug(f"proj_out_all before conv shape: {proj_out_all.shape}")

        # RSD: No need to unsqueeze here. Have batch, channels, height, width
        # proj_out_all = torch.unsqueeze(proj_out_all, 0)
        proj_out_all = torch.permute(proj_out_all, (0, 3, 1, 2))
        proj_out_all = torch.nn.functional.conv2d(
            proj_out_all, filter_2D, padding="same", groups=nPages
        )
        # proj_out_all = torch.squeeze(torch.permute(proj_out_all, (0, 2, 3, 1)), 0)
        proj_out_all = torch.permute(proj_out_all, (0, 2, 3, 1))
        logging.debug(f"proj_out_all after conv shape: {proj_out_all.shape}")

    return proj_out_all


def array_interpolate(
    proj_out_all,
    tomo_obj_all,
    Ax,
    Ay,
    Tx,
    Ty,
    nCols,
    nRows,
):
    """Interpolation with array operations to avoid for loop."""

    device = proj_out_all.device
    logging.debug("Ax shape: {}".format(Ax.shape))

    # RSD: Something weird here. Now the indices are 2D? Will have to debug later. Easiest would probably be to pick fourth index.
    indices = np.argwhere((Ax > 0) & (Ax < nCols) & (Ay > 0) & (Ay < nRows))
    arg_n, arg_0, arg_1, arg_2 = (
        indices[:, 0],
        indices[:, 1],
        indices[:, 2],
        indices[:, 3],
    )

    # RSD: Which to unsqueeze? 1 or 2? Assume 2. HMMM. 2 gives error.
    temp1 = (
        torch.from_numpy(
            Tx[arg_n, arg_0, arg_1, arg_2] * Ty[arg_n, arg_0, arg_1, arg_2]
        ).unsqueeze(1)
    ).to(device)
    temp2 = (
        torch.from_numpy(
            Tx[arg_n, arg_0, arg_1, arg_2] * (1 - Ty[arg_n, arg_0, arg_1, arg_2])
        ).unsqueeze(1)
    ).to(device)
    temp3 = (
        torch.from_numpy(
            (1 - Tx[arg_n, arg_0, arg_1, arg_2]) * Ty[arg_n, arg_0, arg_1, arg_2]
        ).unsqueeze(1)
    ).to(device)
    temp4 = (
        torch.from_numpy(
            (1 - Tx[arg_n, arg_0, arg_1, arg_2]) * (1 - Ty[arg_n, arg_0, arg_1, arg_2])
        ).unsqueeze(1)
    ).to(device)

    logging.debug("indices shape: {}".format(indices.shape))
    logging.debug("tomo_obj_all shape: {}".format(tomo_obj_all.shape))

    # RSD: The issue has been that indexing is occuring in parallel, not accumulative. Hence, a special formula is needed.
    Ax = torch.from_numpy(Ax.astype(np.int64)).to(device)
    Ay = torch.from_numpy(Ay.astype(np.int64)).to(device)
    arg_n = torch.from_numpy(arg_n.astype(np.int64)).to(device)

    proj_out_all = torch.index_put(
        proj_out_all,
        (arg_n, Ay[arg_n, arg_0, arg_1, arg_2], Ax[arg_n, arg_0, arg_1, arg_2]),
        tomo_obj_all[arg_n, arg_0, arg_1, arg_2, :] * temp1,
        accumulate=True,
    )
    proj_out_all = torch.index_put(
        proj_out_all,
        (arg_n, Ay[arg_n, arg_0, arg_1, arg_2] - 1, Ax[arg_n, arg_0, arg_1, arg_2]),
        tomo_obj_all[arg_n, arg_0, arg_1, arg_2, :] * temp2,
        accumulate=True,
    )
    proj_out_all = torch.index_put(
        proj_out_all,
        (arg_n, Ay[arg_n, arg_0, arg_1, arg_2], Ax[arg_n, arg_0, arg_1, arg_2] - 1),
        tomo_obj_all[arg_n, arg_0, arg_1, arg_2, :] * temp3,
        accumulate=True,
    )
    proj_out_all = torch.index_put(
        proj_out_all,
        (arg_n, Ay[arg_n, arg_0, arg_1, arg_2] - 1, Ax[arg_n, arg_0, arg_1, arg_2] - 1),
        tomo_obj_all[arg_n, arg_0, arg_1, arg_2, :] * temp4,
        accumulate=True,
    )

    return proj_out_all


def batch_multiply(A: torch.tensor, B: torch.tensor):
    # Use Einstein sum in different cases
    # Needs improvements for batch, orientation, and more coefficients
    # RSD: Would assume conditions work now.

    # Need possibly one for Rstr * q_obj in batch if matlab does not work.

    if A.shape[-1] == B.shape[-1] and len(B.shape) == 4:  # 3:
        # RSD: Eg sumlm_alm_Ylm
        # RSD: Change condition. Possibly correct now.
        if len(A.shape) == 1:
            A = A.unsqueeze(0).unsqueeze(0)  # RSD: Add two dimensions
        C = torch.einsum("ijm,njkm->nikm", A, B)

    elif len(A.shape) == len(B.shape):  # len(A.shape) > len(B.shape):
        # RSD: No longer the case. Remember to change. If one case per statement, use string.
        # RSD: Eg q_pp

        # C = torch.einsum("ijk,jm->imk", A, B) OLD
        C = torch.einsum("ijk,njm->nimk", A, B)  # RSD: Add one dimension

    elif len(A.shape) < len(B.shape) and len(A.shape) <= 2:
        # RSD: Y_coeff unchanged, num of dims 2 or 1.
        # RSD: Eg Ylm

        if len(A.shape) == 1:
            A = A.unsqueeze(0)  # RSD: Add batch dimension

        C = torch.einsum("ij,njkm->nikm", A, B)
    else:
        raise ValueError(
            "A and B have incompatible shapes. Please implement desired einsum"
        )
    return C


def repmat_cumprod_SH(cos_theta_sh_cut, numOfCoeffs, n_proj):
    """
    Batched version of cumprod_SH
    """
    device = cos_theta_sh_cut.device
    copy_matrix = torch.ones(
        (n_proj, numOfCoeffs, cos_theta_sh_cut.shape[-2], cos_theta_sh_cut.shape[-1]),
        dtype=torch.float64,
    ).to(device)
    # RSD: Add batch dimension.

    for i in range(1, numOfCoeffs):

        copy_matrix[:, i, :, :] = cos_theta_sh_cut ** (2 * i)

    return copy_matrix


def norm_integrand(x, B):

    pi = torch.acos(torch.zeros(1)).item() * 2
    B = B.unsqueeze(0)
    sinx = torch.sin(x)
    while len(sinx.shape) < len(B.shape):
        sinx = sinx.unsqueeze(-1)
    return 2 * pi * sinx * torch.exp(-B * sinx**2)


def torch_simpsons(func, B, a=0, b=torch.acos(torch.zeros(1)).item() * 2, N=500):
    """
    Performs the Simpson's rule integration of a function func
    between a and b with N intervals.
    """
    assert N % 2 == 0, "N must be even"
    h = (b - a) / N
    x = torch.linspace(a, b, N + 1).to(B.device)
    return np.squeeze(
        h
        / 3
        * (
            func(x[0], B)
            + 4 * torch.sum(func(x[1:-1:2], B), dim=0)
            + 2 * torch.sum(func(x[2:-1:2], B), dim=0)
            + func(x[-1], B)
        )
    )


if __name__ == "__main__":

    new_model = 0

    if new_model:
        # RSD: Wrong file.
        workspace = scipy.io.loadmat(
            r"C:\Users\Bruker\OneDrive\Dokumenter\NTNU\XRD_CT\Data sets\Debug Data\workspace_batch_python_orientation.mat"
        )

        theta = workspace["theta_struct"]
        phi = workspace["phi_struct"]
        A = workspace["A"]
        B = workspace["B"]
        unit_q_beamline = workspace["unit_q_beamline"]
        p = workspace["p"]
        X = workspace["X"]
        Y = workspace["Y"]
        Z = workspace["Z"]
        ny = workspace["ny"]
        nx = workspace["nx"]
        nz = workspace["nz"]
        numOfsegments = workspace["numOfsegments"]
        numOfpixels = workspace["numOfpixels"]
        numOfvoxels = workspace["numOfvoxels"]
        find_grad = workspace["find_grad"]

        tic = time.time()

        EXPSIN_main(
            theta,
            phi,
            A,
            B,
            unit_q_beamline,
            p,
            X,
            Y,
            Z,
            ny,
            nx,
            nz,
            numOfsegments,
            numOfpixels,
            numOfvoxels,
            find_grad,
        )

        tac = time.time()
        print("Time elapsed: {}".format(tac - tic))

    else:
        # # Test code
        # # Workspace has to be updated with p.projection_filename.
        # workspace = scipy.io.loadmat(
        #     r"C:\Users\Bruker\OneDrive\Dokumenter\NTNU\XRD_CT\Data sets\Debug Data\workspace_batch_python_orientation.mat"
        # )
        # theta_struct_it = workspace["theta_struct"]
        # phi_struct_it = workspace["phi_struct"]
        # a_temp_it = workspace["a_temp"]
        # ny = workspace["ny"]
        # nx = workspace["nx"]
        # nz = workspace["nz"]
        # numOfsegments = workspace["numOfsegments"]
        # current_projection = (
        #     1  # workspace["projection"][0]  # [2] Load from disk took additional 300ms.
        # )
        # p = workspace["p"]
        # # p[
        # #     "projection_filename"
        # # ] = r"C:\Users\Bruker\OneDrive\Dokumenter\NTNU\XRD_CT\Data sets\Validation_periodic_filter1_3cube_4off_0align_stripped.mat"
        # X = workspace["X"]
        # Y = workspace["Y"]
        # Z = workspace["Z"]
        # numOfpixels = workspace["numOfpixels"]
        # unit_q_beamline = workspace["unit_q_beamline"]
        # Ylm_coef = workspace["Ylm_coef"]
        # find_coefficients = workspace["find_coefficients"]
        # find_orientation = workspace["find_orientation"]
        # numOfCoeffs = workspace["numOfCoeffs"]
        # numOfvoxels = workspace["numOfvoxels"]

        # tic = time.time()

        # SH_main(
        #     theta_struct_it,
        #     phi_struct_it,
        #     a_temp_it,
        #     ny,
        #     nx,
        #     nz,
        #     numOfsegments,
        #     current_projection,
        #     p,
        #     X,
        #     Y,
        #     Z,
        #     numOfpixels,
        #     unit_q_beamline,
        #     Ylm_coef,
        #     find_coefficients,
        #     find_orientation,
        #     numOfCoeffs,
        #     numOfvoxels,
        #     find_grad=True,
        # )
        # tac = time.time()
        # print(f"Time elapsed: {tac - tic}")

        # # Run grad performance test.
        import json

        try:
            with open("performance_results.json") as f:
                performance_results = json.load(f)
        except:
            performance_results = {}

        workspace_names = [
            r"C:\Users\Bruker\OneDrive\Dokumenter\NTNU\XRD_CT\Data sets\Dummy_periodic_filter1_5cube_0off_0alignworkspace.mat",
            r"C:\Users\Bruker\OneDrive\Dokumenter\NTNU\XRD_CT\Data sets\Dummy_periodic_filter1_10cube_0off_0alignworkspace.mat",
            r"C:\Users\Bruker\OneDrive\Dokumenter\NTNU\XRD_CT\Data sets\Dummy_periodic_filter1_11cube_4off_0alignworkspace.mat",
            r"C:\Users\Bruker\OneDrive\Dokumenter\NTNU\XRD_CT\Data sets\Dummy_periodic_filter1_20cube_0off_0alignworkspace.mat",
            r"C:\Users\Bruker\OneDrive\Dokumenter\NTNU\XRD_CT\Data sets\Dummy_periodic_filter1_17cube_8off_0alignworkspace.mat",
            r"C:\Users\Bruker\OneDrive\Dokumenter\NTNU\XRD_CT\Data sets\Dummy_periodic_filter1_30cube_0off_0alignworkspace.mat",
            r"C:\Users\Bruker\OneDrive\Dokumenter\NTNU\XRD_CT\Data sets\Dummy_periodic_filter1_27cube_8off_0alignworkspace.mat",
            r"C:\Users\Bruker\OneDrive\Dokumenter\NTNU\XRD_CT\Data sets\Dummy_periodic_filter1_37cube_8off_0alignworkspace.mat",
        ]
        workspace_names = [
            r"/home/rubensd/Documents/MATLAB/Data sets/Dummy_periodic_filter1_5cube_0off_0alignworkspace.mat",
            r"/home/rubensd/Documents/MATLAB/Data sets/Dummy_periodic_filter1_10cube_0off_0alignworkspace.mat",
            r"/home/rubensd/Documents/MATLAB/Data sets/Dummy_periodic_filter1_11cube_4off_0alignworkspace.mat",
            r"/home/rubensd/Documents/MATLAB/Data sets/Dummy_periodic_filter1_20cube_0off_0alignworkspace.mat",
            r"/home/rubensd/Documents/MATLAB/Data sets/Dummy_periodic_filter1_17cube_8off_0alignworkspace.mat",
            r"/home/rubensd/Documents/MATLAB/Data sets/Dummy_periodic_filter1_30cube_0off_0alignworkspace.mat",
            r"/home/rubensd/Documents/MATLAB/Data sets/Dummy_periodic_filter1_27cube_8off_0alignworkspace.mat",
            r"/home/rubensd/Documents/MATLAB/Data sets/Dummy_periodic_filter1_37cube_8off_0alignworkspace.mat",
        ]
        voxels = np.array(
            [5**3, 10**3, 15**3, 20**3, 25**3, 30**3, 35**3, 45**3]
        )
        gradient_compute_times = np.zeros(len(workspace_names))

        simulation = "GPU"
        performance_results[simulation] = {}

        performance_results[simulation]["workspace_names"] = workspace_names
        performance_results[simulation]["voxels"] = voxels

        for i, path in enumerate(workspace_names):
            workspace = scipy.io.loadmat(path)
            theta_struct_it = workspace["theta_struct"]
            phi_struct_it = workspace["phi_struct"]
            a_temp_it = workspace["a_temp"]
            ny = workspace["ny"]
            nx = workspace["nx"]
            nz = workspace["nz"]
            numOfsegments = workspace["numOfsegments"]
            current_projection = path[: -(len("workspace.mat"))] + ".mat"
            print(current_projection)
            p = workspace["p"]
            X = workspace["X"]
            Y = workspace["Y"]
            Z = workspace["Z"]
            numOfpixels = workspace["numOfpixels"]
            unit_q_beamline = workspace["unit_q_beamline"]
            Ylm_coef = workspace["Ylm_coef"]
            find_coefficients = 1  # workspace["find_coefficients"]
            find_orientation = 1  # workspace["find_orientation"]
            numOfCoeffs = workspace["numOfCoeffs"]
            numOfvoxels = workspace["numOfvoxels"]

            tic = time.time()

            SH_main(
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
            gradient_compute_times[i] = tac - tic

        performance_results[simulation][
            "gradient_compute_times"
        ] = gradient_compute_times

        try:
            with open("performance_results.json", "w") as f:
                json.dump(performance_results, f)
        except:
            print("Failed to save performance results.")
            np.save(f"{simulation}_performance_results.npy", performance_results)
            np.save(f"{simulation}_gradient_compute_times.npy", gradient_compute_times)
        finally:
            print("Performance results saved.")
