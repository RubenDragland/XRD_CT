

function [E, grad_a, grad_theta_struct, grad_phi_struct, aux_diff_poisson, proj_out_all ] = SAXS_AD_all_forward_backward(theta_struct, phi_struct, a_temp, ny, nx, nz, numOfsegments, projection, p, X, Y, Z, numOfpixels, unit_q_beamline, Ylm_coef)
    
    E = 0; %RSD: Init error. 
    grad_a = zeros(ny, nx, nz, numOfCoeffs); %RSD: Init grads
    grad_theta_struct = zeros(ny, nx, nz);
    grad_phi_struct   = zeros(ny, nx, nz);

    parfor ii = 1: length(projection)

        current_projection = projection(ii);

        if find_coefficients
            a_temp = dlarray(a_temp);
        end

        theta_struct = dlarray(theta_struct);
        phi_struct = dlarray(phi_struct);

        if find_coefficients %RSD: Consider to remove this if/else. 
            error_norm, AD_grad_coeffs, AD_grad_theta, AD_grad_phi, aux_diff_poisson, proj_out_all = dlfeval(@SAXS_AD_all_cost_function, theta_struct, phi_struct, a_temp, ny, nx, nz, numOfsegments, current_projection, p, X, Y, Z, numOfpixels, unit_q_beamline, Ylm_coef, find_coefficients );
            
            AD_grad_coeff = extractdata(AD_grad_coeff); % RSD: No longer need for dlarrays
            grad_a = grad_a + reshape(AD_grad_coeff, ny, nx, nz, numOfCoeffs);
        else
            error_norm, AD_grad_theta, AD_grad_phi, aux_diff_poisson, proj_out_all = dlfeval(@SAXS_AD_all_cost_function, theta_struct, phi_struct, a_temp, ny, nx, nz, numOfsegments, current_projection, p, X, Y, Z, numOfpixels, unit_q_beamline, Ylm_coef, find_coefficients);
        end

        error_norm = extractdata(error_norm);
        E = E + error_norm;

        aux_diff_poisson = extractdata(aux_diff_poisson);
        proj_out_all = extractdata(proj_out_all);

        AD_grad_theta = extractdata(AD_grad_theta);
        grad_theta_struct = grad_theta_struct + reshape(AD_grad_theta, ny, nx, nz); 
        
        AD_grad_phi = extractdata(AD_grad_phi);    
        grad_phi_struct = grad_phi_struct + reshape(AD_grad_phi, ny, nx, nz); 

        return_synth_proj = 0; %RSD: Manual debug param.

        if return_synth_proj
            proj_out(ii).projection = proj_out_all;
            proj_out(ii).errorplot = aux_diff_poisson; % Spatially resolved error
            proj_out(ii).error = error_norm;           % error metric with poisson noise
            proj_out(ii).rotx = projection(ii).rot_x;
            proj_out(ii).roty = projection(ii).rot_y;
        end

    end

end


function [ error_norm, ad_grad_coeffs, ad_grad_theta, ad_grad_phi, aux_diff_poisson, proj_out_all] = SAXS_AD_all_cost_function(theta_struct, phi_struct, a_temp, ny, nx, nz, numOfsegments, current_projection, p, X, Y, Z, numOfpixels, unit_q_beamline, Ylm_coef, find_coefficients )
    
    
    data = double(current_projection.data) ; %RSD: Change relative to the other cost function.
    Rot_exp_now = double(current_projection.Rot_exp) ; 
    unit_q_object = Rot_exp_now' * unit_q_beamline;

    sin_theta_struct = reshape(sin(theta_struct), 1, 1, numOfvoxels);
    cos_theta_struct = reshape(cos(theta_struct), 1, 1, numOfvoxels);
    sin_phi_struct = reshape(sin(phi_struct), 1, 1, numOfvoxels);
    cos_phi_struct = reshape(cos(phi_struct), 1, 1, numOfvoxels);
    
    Rot_str = [ ...
        cos_theta_struct.*cos_phi_struct, cos_theta_struct.*sin_phi_struct, -sin_theta_struct; ...
        -sin_phi_struct                 , cos_phi_struct                  , zeros_struct     ; ...
        sin_theta_struct.*cos_phi_struct, sin_theta_struct.*sin_phi_struct, cos_theta_struct];

    
    q_pp = pagemtimes(Rot_str, unit_q_object); %bsxpagemult(double(Rot_str), double(unit_q_object)); %RSD: Change to pagemtimes. 

    cos_theta_sh_cut = q_pp(3, :, :); 

    if all(mod(l,2) == 0) && all(m == 0) % Make sure that all l is even, and that all m = 0
        block_cos_theta_powers = double([ones_struct; cumprod(repmat(cos_theta_sh_cut.^2, numOfCoeffs-1, 1), 1)]); %RSD: Check functionality of cumprod, repmat etc.  % RSD: ONES AND CUMULATIVE PRODUCT OF CREATED NDARRAY.
    elseif ~all(mod(l,2) == 0)
        error('Not all orders l of SH are even, odd orders are currently not supported')
    elseif ~all(m == 0)
        error('Not all orders m of SH are zero, only m = 0 are currently supported')
    end

    Ylm = pagemtimes(Ylm_coef, block_cos_theta_powers); %bsxpagemult(double(Ylm_coef), block_cos_theta_powers); %RSD: change to pagemtimes.

    sumlm_alm_Ylm = pagemtimes(a_temp, Ylm); %bsxpagemult(double(a_temp), Ylm); %RSD: Change to pagemtimes

    data_synt_vol = permute(abs(sumlm_alm_Ylm.^2), [3, 2, 1]); % modeled intensities
    data_synt_vol = reshape(data_synt_vol, ny, nx, nz, numOfsegments);
    
    %%%% Projection of the data_synt_vol onto detector plane, in direction
    %%%% of Rot_exp
    
    % output of projection should have the size of the measured projection
    %put in alignment correction values from registration
    xout = (1:size(data,2)) - ceil(size(data,2)/2) + current_projection.dx;     % RSD: DATA IS JUST PROJECTION DATA
    yout = (1:size(data,1)) - ceil(size(data,1)/2) + current_projection.dy;
    
    % I [y x numOfsegments] Eq(3) of Liebi et al
    [proj_out_all, xout, yout] = arb_projection(data_synt_vol, X, Y, Z, Rot_exp_now, p, xout, yout);   % RSD: GET ESTIMATED INTENSITIES. %RSD: Should work without tweaking. 
    data_synt_vol = []; % free up memory

    aux_diff_poisson = (sqrt(proj_out_all) - sqrt(data)) .* current_projection.window_mask; %projection(ii).window_mask;
    error_norm = 2*sum(sum(sum(aux_diff_poisson.^2))) / numOfpixels; % error metric with poisson noise

    if find_coefficients
        try
            ad_grad_coeffs = dlgradient(error_norm, a_temp);
        catch
            ad_grad_coeffs = dlarray(zeros(size(a_temp)));
        end
    else
        try
            ad_grad_theta = dlgradient(error_norm, theta_struct);
            ad_grad_phi = dlgradient(error_norm, phi_struct);
        catch
            ad_grad_theta = dlarray(zeros(size(theta_struct)));
            ad_grad_phi = dlarray(zeros(size(phi_struct)));
        end
    end

end