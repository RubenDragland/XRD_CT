
% RSD: The following functions are based on the codeblock in the error metric script, with AD implementation as an addition. 


function [E, AD_grad_coeff] = SAXS_AD_optim_fb_coefficients(a_temp, theta_struct, phi_struct, ny, nx, nz, numOfsegments, projection, p, X, Y, Z, numOfpixels, unit_q_beamline, Ylm_coef, find_coefficients, numOfCoeffs, numOfvoxels)

    sin_theta_struct = reshape(sin(theta_struct), 1, 1, numOfvoxels);
    cos_theta_struct = reshape(cos(theta_struct), 1, 1, numOfvoxels);
    sin_phi_struct = reshape(sin(phi_struct), 1, 1, numOfvoxels);
    cos_phi_struct = reshape(cos(phi_struct), 1, 1, numOfvoxels);

    %%% Define the rotation matrices
    % Rot_str is the combination of 1. Rotation around y with theta_struct
    % and 2. Rotation around z with theta_struct
    % this rotation is needed from the spherical harmonics coordinate system (main
    % orientation axis along z-axis) to the object coordinate system
    Rot_str = [ ...
        cos_theta_struct.*cos_phi_struct, cos_theta_struct.*sin_phi_struct, -sin_theta_struct; ...
        -sin_phi_struct                 , cos_phi_struct                  , zeros_struct     ; ...
        sin_theta_struct.*cos_phi_struct, sin_theta_struct.*sin_phi_struct, cos_theta_struct];


    %calculate for all projections
    parfor ii = 1:length(projection) %use parallel processing for the loop over all projections
    %for ii = 1:length(projection) % RSD: Debug for loop  
        
        data = double(projection(ii).data);
        
        Rot_exp_now = double(projection(ii).Rot_exp);       % Rotation matrix of projection (3x3), R_exp

        % Coordinates of the measured intensities for each sector, in reciprocal space and in object coordinates.  
        % q' = R_exp * q (unit vector)
        % Size 3x8 [xyz  numOfsegments]
        unit_q_object = Rot_exp_now' * unit_q_beamline;      % RSD: Accounting for the rotation of the projection to the beamline

        % Coordinates xyz of the measured intensities, for each sector and each
        % voxel in reciprocal space and in coordinates centered on the
        % spherical harmonic local orientation vector.
        %   q'' = Rot_str*q'      This q'' coordinates are not defined in the paper
        % Liebi et al 2018 as such, it is another coordinate system aligned to
        % the main axes of the SH functions, basically aligned to u_str
        %   [3 numOfsegments numOfvoxels] = [3 3 numOfvoxels]*[3 numOfsegments]
        q_pp = bsxpagemult(double(Rot_str), double(unit_q_object)); % RSD: Matrix multiplication. Using rotation matrix. q in terms of SH
        
        % The third component of this q'' vector corresponds to the cosine of
        % the theta angle that defines the location of the measurement in terms
        % of spherical harmonics
        %   [1xnumOfsegments numOfvoxels]
        cos_theta_sh_cut = q_pp(3, :, :); 
        

        block_cos_theta_powers = double([ones_struct; cumprod(repmat(cos_theta_sh_cut.^2, numOfCoeffs-1, 1), 1)]);    % RSD: ONES AND CUMULATIVE PRODUCT OF CREATED NDARRAY.

        % After multiplication with the coefficients of the Legendre polynomials Ylm becomes the properly normalized SH functions
        %   [numOforder numOfsegments numOfvoxels] = [numOforder numOforder]*[numOforder numOfsegments numOfvoxels]    .
        Ylm = bsxpagemult(double(Ylm_coef), block_cos_theta_powers);  % RSD: RETRIEVE THE SH FUNCTION FROM MATRIX PRODUCT WITH ARRAY OF POWERS OF COSINE. YLM-COEFFICIENTS MATCHED WITH ITS "COSINE-POLYNOMIAL"

        
        if find_grad || p.GPU % RSD: Need to keep Matlab if finding gradient.
            sumlm_alm_Ylm = pagemtimes(a_temp, Ylm); % RSD: MATLAB option for trace of dlarray?
        else
            sumlm_alm_Ylm = bsxpagemult(double(a_temp), Ylm );   % RSD: INNER PART OF EQ(3) OR EQ(1) LIEBI ET AL 2018
        end
    

        data_synt_vol = permute(abs(sumlm_alm_Ylm.^2), [3, 2, 1]); 
        data_synt_vol = reshape(data_synt_vol, ny, nx, nz, numOfsegments);
    
    
        % output of projection should have the size of the measured projection
        %put in alignment correction values from registration
        xout = (1:size(data,2)) - ceil(size(data,2)/2) + current_projection.dx; %projection(ii).dx;     % RSD: DATA IS JUST PROJECTION DATA
        yout = (1:size(data,1)) - ceil(size(data,1)/2) + current_projection.dy; %projection(ii).dy;     % RSD: Parfor overhead communication issue lead to this change.
    
        % I [y x numOfsegments] Eq(3) of Liebi et al
        [proj_out_all, xout, yout] = arb_projection(data_synt_vol, X, Y, Z, Rot_exp_now, p, xout, yout, find_grad);   % RSD: GET ESTIMATED INTENSITIES. May impose an issue.
        data_synt_vol = []; % free up memory
    
        %%%calculating the error (difference between estimated intensity and measured data), taking poisson noise into account
        %auxillary functions, used in E and the gradient
        aux_diff_poisson = (sqrt(proj_out_all) - sqrt(data)) .* current_projection.window_mask; %projection(ii).window_mask;
        error_norm = 2*sum(sum(sum(aux_diff_poisson.^2))) / numOfpixels; % error metric with poisson noise

        E = E + error_norm;

    end
    
    if find_grad
        %try
        %    AD_grad_coeff = dlgradient(E, a_temp);
        %catch
        %    AD_grad_coeff = dlarray( zeros(size(a_temp)) ); % RSD: Assuming the gradients are 0 if none elements of original dlarray are transferred to the projection.
        %end
        AD_grad_coeff = dlgradient(E, a_temp); 
    else
        AD_grad_coeff = [];   %RSD: Return empty if no gradient is to be calculated.
    end
end
