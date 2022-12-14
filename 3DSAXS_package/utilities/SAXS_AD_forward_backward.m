
% RSD: The following functions are based on the codeblock in the error metric script, with AD implementation as an addition. 


function [error_norm, ad_grad ] = SAXS_AD_forward_backward(func, a_temp, Ylm, ny, nx, nz, numOfsegments, data, current_projection, Rot_exp_now, p, find_grad, X, Y, Z, numOfpixels ) % RSD: EXPAND TO OTHER UTILITIES LATER.
    if find_grad
        a_temp = dlarray(a_temp); % RSD: Also fix function call. %RSD: gpuArray( dlarray(a_temp ) ); Need to get it working on server. 
        if p.GPU
            a_temp = gpuArray(a_temp); %RSD: Have this for only find grad or both?
        end
        [error_norm, ad_grad] = dlfeval(@SAXS_AD_cost_function, a_temp, Ylm, ny, nx, nz, numOfsegments, data, current_projection, Rot_exp_now, p, find_grad, X, Y, Z, numOfpixels ); %, Ylm, ny, nx, nz, numOfsegments, data, current_projection, Rot_exp_now, p, find_grad, X, Y, Z, numOfpixels);
        if p.GPU
            ad_grad = gather(ad_grad);
            error_norm = gather(error_norm);
        end
    else
        
        [error_norm, ad_grad] = SAXS_AD_cost_function( a_temp, Ylm, ny, nx, nz, numOfsegments, data, current_projection, Rot_exp_now, p, find_grad, X, Y, Z, numOfpixels ); %, Ylm, ny, nx, nz, numOfsegments, data, current_projection, Rot_exp_now, p, find_grad, X, Y, Z, numOfpixels);
        
    end
    % RSD: Save memory by using nested functions. Fewer input params. Not allowed for parallel loop. 
end

function [error_norm, ad_grad ] = SAXS_AD_cost_function(a_temp, Ylm, ny, nx, nz, numOfsegments, data, current_projection, Rot_exp_now, p, find_grad, X, Y, Z, numOfpixels ) % Ylm, ny, nx, nz, numOfsegments, data, current_projection, Rot_exp_now, p, find_grad, X, Y, Z, numOfpixels)

    %sumlm_alm_Ylm = bsxpagemult(a_temp, Ylm);
    if find_grad || p.GPU % RSD: Need to keep Matlab if finding gradient.
        sumlm_alm_Ylm = pagemtimes(a_temp, Ylm); % RSD: MATLAB option for trace of dlarray?
    else
        sumlm_alm_Ylm = bsxpagemult(double(a_temp), Ylm );   % RSD: INNER PART OF EQ(3) OR EQ(1) LIEBI ET AL 2018
    end


    %RSD: BELIEVE THE ORDERING BECOMES:
    %{1 NUM_SEGMENTS NUM_VOXELS] = [1 NUM_VOXELS NUM_ORDER/COEFF] X [NUM_ORDER NUM_SEGMENTS NUM_VOXELS] ?

    % R [y x z numOfsegments]   Eq. (1) of Liebi et al 2018, here organized
    % as [1 numOfsegments numOfvoxels]
    data_synt_vol = permute(abs(sumlm_alm_Ylm.^2), [3, 2, 1]); % modeled intensities %RSD: Check expression.
    data_synt_vol = reshape(data_synt_vol, ny, nx, nz, numOfsegments);

    %%%% Projection of the data_synt_vol onto detector plane, in direction
    %%%% of Rot_exp

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

    if find_grad
        try
            ad_grad = dlgradient(error_norm, a_temp);
        catch
            ad_grad = dlarray( zeros(size(a_temp)) ); % RSD: Assuming the gradients are 0 if none elements of original dlarray are transferred to the projection.
        end
    else
        ad_grad = [];   %RSD: Return empty if no gradient is to be calculated.
    end
    % RSD: Ignore so we may locate crash.
    % RSD: BACKPROJECTION OF THE GRADIENT. However, should one possibly use aux_diff_poisson as function instead. The number is really not that useful.
    % RSD: HOWEVER, IN THEORY, THIS IS NO BIG DEAL. THE QUESTION IS THE OUTPUT OF DLGRADIENT. 
    end

