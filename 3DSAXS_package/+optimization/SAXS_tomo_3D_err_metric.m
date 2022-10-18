function [E, grad, proj_out, Ereg] = SAXS_tomo_3D_err_metric(opt_inputs, p, s, projection)
% ** opt_inputs     Vector with input variables under optimization, set to
%                   empty to skip optimization and just compute error metric 
%                   or provide an output projection for data evaluation.
% 
% returns:
% ++ E          Error metric
% ++ grad       Gradient of E with respect to optimization variables
% ++ proj_out   Synthetic projections from the 3D data, also includes the
%               2D error map and the data error per projection
% ++ Ereg       Orientation regularization error metric

%*-------------------------------------------------------------------------------------*
%|                                                                                     |
%|  Except where otherwise noted, this work is licensed under a                        |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0                          |
%|  International (CC BY-NC-SA 4.0) license.                                           |
%|                                                                                     |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)                  |
%|                                                                                     |
%|      Author: CXS group, PSI                                                         |
%*------------------------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form and additionally you should cite:
%     M. Liebi, M. Georgiadis, A. Menzel, P. Schneider, J. Kohlbrecher, 
%     O. Bunk, and M. Guizar-Sicairos, “Nanostructure surveys of 
%     macroscopic specimens by small-angle scattering tensor tomography,”
%     Nature 527, 349-352 (2015).   (doi:10.1038/nature16056)
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of coefficients
numOfCoeffs = numel([s.a.l]);

if isempty(opt_inputs)
    skip_optimization = true;
else
    skip_optimization = false;
end

if numOfCoeffs > 4
    error('add more coefficients in error metric');
end

% start a parallel pool (if none exists yet) with a number of workers equal to the number
% of cores on the current machine 
if isempty(gcp('nocreate'))
    parpool('local', feature('numcores'));
end

% parameters for the optimization
return_synth_proj = (nargout>=3);   % Return projections from synthetic data in proj_out
return_Ereg = (nargout>=4);         % Return value of orientation regularization
if skip_optimization
    warning('Skipping optimization')
    find_coefficients = false;
    find_orientation  = false;
    grad = [];
    find_grad = false;
    p.opt_coeff(:) = false;
else
    find_coefficients = any(p.opt_coeff);       % Optimize over SH coefficients: a (true or false);
    find_orientation = p.find_orientation;
    find_grad = (nargout >= 2);
    if (p.regularization_angle) && (~find_orientation)
        warning('p.regularization_angle was true but optimization of orientation (p.find_orientation) is false. Regularization of angle will be turned off as it does not make sense if the orienation is not optimized')
        p.regularization_angle = 0; % Notice this still allows above to compute the regularization skip_optimization
    end
end

phi_det = p.phi_det;                            % Measurement azimuthal angles (detector)
theta_det = ones(size(phi_det)).*projection(1).integ.theta_det;     
% Polar angle of SH in beamline coordinates, e.g. = pi/2 + Ewald sphere correction
l = [s.a.l];                                        % Polar order of the spherical harmonics
m = [s.a.m];

nx = p.nx;                            
ny = p.ny;
nz = p.nz;
numOfsegments = p.numsegments;                      % number of segements the saxs data is averaged (default is 16)
numOfvoxels = p.numOfvoxels;            % dimension of tomogram
mask3D = s.mask3D;
numOfpixels = numel(projection)*nx*ny;
numOfordersopt = sum(logical(p.opt_coeff)); % Number of coefficients that are being optimized

%opt_inputs is the optimization vector, containing theta_struct, phi_struct and a, the coefficents of the spherical harmonics
if find_orientation
    theta_struct = opt_inputs(1:numOfvoxels);
    phi_struct = opt_inputs(numOfvoxels+1:2*numOfvoxels);
    theta_struct = reshape(theta_struct,[p.ny,p.nx,p.nz]);
    phi_struct   = reshape(phi_struct  ,[p.ny,p.nx,p.nz]); 
    opt_inputs = opt_inputs(2*numOfvoxels+1:end);
    
    %avoid wrapping of angles
    if p.avoid_wrapping
        phi_struct = mod(phi_struct, 2*pi);
        theta_struct = mod(theta_struct, 2*pi);
    end
else
    theta_struct = s.theta.data;
    phi_struct = s.phi.data;
end

a = [];
for ii = 1:numel(p.opt_coeff)
    if p.opt_coeff(ii)
        a = [a opt_inputs(1:numOfvoxels)];
        opt_inputs = opt_inputs(numOfvoxels+1:end);
    else
        a = [a s.a(ii).data(:).'];
    end
end
a = reshape(a, numOfvoxels, numOfCoeffs);

%initalize arrays for the gradients
grad_a = zeros(ny, nx, nz, numOfCoeffs);

grad_theta_struct = zeros(ny, nx, nz);
grad_phi_struct   = zeros(ny, nx, nz);

sin_theta_struct = reshape(sin(theta_struct), 1, 1, numOfvoxels);
cos_theta_struct = reshape(cos(theta_struct), 1, 1, numOfvoxels);
sin_phi_struct = reshape(sin(phi_struct), 1, 1, numOfvoxels);
cos_phi_struct = reshape(cos(phi_struct), 1, 1, numOfvoxels);

zeros_struct = zeros(1, 1, numOfvoxels);

%%% Define the rotation matrices
% Rot_str is the combination of 1. Rotation around y with theta_struct
% and 2. Rotation around z with theta_struct
% this rotation is needed from the spherical harmonics coordinate system (main
% orientation axis along z-axis) to the object coordinate system
Rot_str = [ ...
    cos_theta_struct.*cos_phi_struct, cos_theta_struct.*sin_phi_struct, -sin_theta_struct; ...
    -sin_phi_struct                 , cos_phi_struct                  , zeros_struct     ; ...
    sin_theta_struct.*cos_phi_struct, sin_theta_struct.*sin_phi_struct, cos_theta_struct];

% the derivative of the rotation matrix with respect to theta_struct and phi_struct
Rot_str_diff_theta = [ ...
    -sin_theta_struct.*cos_phi_struct, -sin_theta_struct.*sin_phi_struct, -cos_theta_struct ; ...
    zeros_struct                     , zeros_struct                     , zeros_struct      ; ...
    cos_theta_struct.*cos_phi_struct , cos_theta_struct.*sin_phi_struct , -sin_theta_struct];

Rot_str_diff_phi = [ ...
    -cos_theta_struct.*sin_phi_struct, cos_theta_struct.*cos_phi_struct, zeros_struct ; ...
    -cos_phi_struct                  , -sin_phi_struct                 , zeros_struct ; ...
    -sin_theta_struct.*sin_phi_struct, sin_theta_struct.*cos_phi_struct, zeros_struct];

% Determine Ylm and Ylm+1 coefficients
Ylm_coef        = optimization.get_Ylm_coeff(l,m);

% Ylm_coef = ...
%     [sqrt(1/pi)*[  1/2,      0,       0,      0]; ...
%     sqrt(5/pi) *[ -1/4,    3/4,       0,      0]; ...
%     sqrt(1/pi) *[ 9/16,  -45/8,  105/16,      0]; ...
%     sqrt(13/pi)*[-5/32, 105/32, -315/32, 231/32]];
% Ylm_coef = Ylm_coef(1:order, 1:order);

% YlmPLUS1_coef = ...
%     [                [0,     0,     0,      0]; ...
%     sqrt(15/(2*pi)) *[0,  -1/2,     0,      0]; ...
%     sqrt(5/pi)      *[0,   9/8, -21/8,      0]; ...
%     sqrt(273/(2*pi))*[0, -5/16,  15/8, -33/16]];
% YlmPLUS1_coef = YlmPLUS1_coef(1:order, 1:order);

a_temp = reshape(a.', 1, numOfCoeffs, []);


a_temp1 = [];
a_temp2 = [];
if find_orientation
    a_temp1 = a .* m * Ylm_coef;
    a_temp1 = reshape(a_temp1.', 1, numOfCoeffs, []);
    
    YlmPLUS1_coef   = optimization.get_Ylm_coeff(l,m+1);
    
    a_temp2 = a .* sqrt((l-m).*(l+m+1)) * YlmPLUS1_coef;
    a_temp2 = reshape(a_temp2.', 1, numOfCoeffs, []);
end

ones_struct = ones(1, numOfsegments, numOfvoxels);

%RSD: SO FAR SO GOOD, MAY PROPABLY BE ABLE TO KEEP EVERYTHING UP TO THIS POINT. 
% optimization
E = 0;


% q (unit vector), 3D coordinates of the measurements in reciprocal space and in beamline coordinates. 
% Size 3x8 [xyz x numOfsegments]. The third dimension is actually only needed to account for Ewald sphere curvature
% correction which is included in theta_det
unit_q_beamline = [sin(theta_det).*cos(phi_det); sin(theta_det).*sin(phi_det);cos(theta_det)]; 

N = [ny nx nz]; % Size of tomogram
x = (1:N(2)) - ceil(N(2)/2);
y = (1:N(1)) - ceil(N(1)/2);
z = (1:N(3)) - ceil(N(3)/2);
[X, Y, Z] = meshgrid(x, y, z); % (x,y,z) (2,1,3)
% If x is defined as the 3rd dimension N(3), then it should appear in the
% third position of the meshgrid command. If x is defined as the first
% dimension N(1), then it appears in the second position of the meshgrid,
% because meshgrid is (x,y,z) and size is (y,x,z)


%calculate for all projections
parfor ii = 1:length(projection) %use parallel processing for the loop over all projections
%for ii = 1:length(projection)
    
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
    
    % Array of powers of cos(theta_sh)^2, important to evaluate the spherical
    % harmonic functions.
    % https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
    %   [numOforder numOfsegments numOfvoxels]
    %   along the first index we have   cos(theta_sh).^(2*[0:order-1])
    % MGS - note that because the cosines are in a particular order and
    % only powers of 2 are included this implies a specific arrangement of
    % coefficients and also excludes m ~= 0
    if all(mod(l,2) == 0) && all(m == 0) % Make sure that all l is even, and that all m = 0
        block_cos_theta_powers = double([ones_struct; cumprod(repmat(cos_theta_sh_cut.^2, numOfCoeffs-1, 1), 1)]);    % RSD: ONES AND CUMULATIVE PRODUCT OF CREATED NDARRAY.
        % As a reminder, when generalizing to any scrambling of coefficient
        % order (e.g. s.a([3 2 1 4]) I also mistakenly changed the order of
        % the powers of cos(theta), this was a mistake because in
        % Ylm_coeff = get_Ylm_coeffs, the output columns are not scrambled
        % but follow the ascending order of coefficients. This relationship
        % between the second index in Ylm_coeff and the power of cos(theta)
        % used needs to be preserved in order for the calculation to be
        % correct. As a reminder, below I put the code that I mistakenly
        % wrote when trying to generalize the coefficient scrambling
        %         block_cos_theta_powers = double(repmat(cos_theta_sh_cut.^2,order,1));
        %         for jj = 1:order
        %              block_cos_theta_powers(jj,:,:) = block_cos_theta_powers(jj,:,:).^(l(jj)/2);
        %         end
    elseif ~all(mod(l,2) == 0)
        error('Not all orders l of SH are even, odd orders are currently not supported')
    elseif ~all(m == 0)
        error('Not all orders m of SH are zero, only m = 0 are currently supported')
    end
    % After multiplication with the coefficients of the Legendre polynomials Ylm becomes the properly normalized SH functions
    %   [numOforder numOfsegments numOfvoxels] = [numOforder numOforder]*[numOforder numOfsegments numOfvoxels]  .
    Ylm = bsxpagemult(double(Ylm_coef), block_cos_theta_powers);  % RSD: RETRIEVE THE SH FUNCTION FROM MATRIX PRODUCT WITH ARRAY OF POWERS OF COSINE. YLM-COEFFICIENTS MATCHED WITH ITS "COSINE-POLYNOMIAL"
    
    % sum_lm(a_lm*Y_lm)
    %   [1 numOfsegments numOfvoxels] = [1 numOforders numOfvoxels] * [numOforder numOfsegments numOfvoxels]
    
    %RSD: START OF ERROR ESTIMATION. A_TEMP HAS TO BE INCLUDED IN
    %COMPUTATIONAL GRAPH. ONE IDEA IS TO IMPLEMENT THIS PART WITH
    %COEFFICIENTS FIRST, AND THEN FIGURE OUT THE ORIENTATION PART
    %MOREOVER, BSXPAGEMULT IS PROPABLY FAST (SHOULD BE CHECKED) THUS NICE
    %TO KEEP IT.
    
    sumlm_alm_Ylm = bsxpagemult(double(a_temp), Ylm);   % RSD: INNER PART OF EQ(3) OR EQ(1) LIEBI ET AL 2018
    %RSD: BELIEVE THE ORDERING BECOMES:
    %{1 NUM_SEGMENTS NUM_VOXELS] = [1 NUM_VOXELS NUM_ORDER/COEFF] X [NUM_ORDER NUM_SEGMENTS NUM_VOXELS] ?
    
    % R [y x z numOfsegments]   Eq. (1) of Liebi et al 2018, here organized
    % as [1 numOfsegments numOfvoxels]
    data_synt_vol = permute(abs(sumlm_alm_Ylm.^2), [3, 2, 1]); % modeled intensities
    data_synt_vol = reshape(data_synt_vol, ny, nx, nz, numOfsegments);
    
    %%%% Projection of the data_synt_vol onto detector plane, in direction
    %%%% of Rot_exp
    
    % output of projection should have the size of the measured projection
    %put in alignment correction values from registration
    xout = (1:size(data,2)) - ceil(size(data,2)/2) + projection(ii).dx;     % RSD: DATA IS JUST PROJECTION DATA
    yout = (1:size(data,1)) - ceil(size(data,1)/2) + projection(ii).dy;
    
    % I [y x numOfsegments] Eq(3) of Liebi et al
    [proj_out_all, xout, yout] = arb_projection(data_synt_vol, X, Y, Z, Rot_exp_now, p, xout, yout);   % RSD: GET ESTIMATED INTENSITIES.
    data_synt_vol = []; % free up memory
    
    %%%calculating the error (difference between estimated intensity and measured data), taking poisson noise into account
    %auxillary functions, used in E and the gradient
    aux_diff_poisson = (sqrt(proj_out_all) - sqrt(data)) .* projection(ii).window_mask;
    error_norm = 2*sum(sum(sum(aux_diff_poisson.^2))) / numOfpixels; % error metric with poisson noise
    E = E + error_norm;
    
    
    %RSD: ERROR IS FOUND. SOME PARTS HAVE TO BE CHANGED AND WRITTEN IN
    %PYTHON IN ORDER TO UTILISE PYTORCH AUTOGRAD.
    
    % When it is requested to have output projection and/or error, useful
    % for debugging
    if return_synth_proj
        proj_out(ii).projection = proj_out_all;
        proj_out(ii).errorplot = aux_diff_poisson; % Spatially resolved error
        proj_out(ii).error = error_norm;           % error metric with poisson noise
        proj_out(ii).rotx = projection(ii).rot_x;
        proj_out(ii).roty = projection(ii).rot_y;
    end
    %%%the gradients RSD: HERE THE GRADIENTS BEGIN.
    if find_grad
        aux_grad_poisson =  aux_diff_poisson ./ (sqrt(proj_out_all) + 1e-10) / numOfpixels; 
        aux_grad_poisson_vol = arb_back_projection(aux_grad_poisson, xout, yout, X, Y, Z, Rot_exp_now, p); 
        
        if find_coefficients
            % = conj( sum(a_lm*Y_lm) )*Ylm
            % [y x z numOfsegments numOforders] = repmat([1 numOfsegments numOfvoxels]) * [numOforders numOfsegments numOf voxels]
            Ymn_aux_vol = permute(real(repmat(conj(sumlm_alm_Ylm), numOfCoeffs, 1, 1) .* Ylm), [3, 2, 1]);
            Ymn_aux_vol = reshape(Ymn_aux_vol, ny, nx, nz, numOfsegments, numOfCoeffs);
            grad_a = grad_a + 4*squeeze(sum(aux_grad_poisson_vol.*Ymn_aux_vol, 4));
            Ymn_aux_vol = []; % free up memory
        end
        Ylm = []; %#ok<*NASGU> % free up memory
        
        
        
        if find_orientation
            % dq''/dthetaop = d(R_str*q')/dthetaop
            % [3 numOfsegments numOfvoxels] = [3 3 numOfvoxels]*[3 numofsegments] 
            % MGS - note, this is not a good name, the output is not a
            % rotation but it is the derivative of measurement coordinates
            % in SH space 
            dqpp_dthetaop = bsxpagemult(Rot_str_diff_theta, unit_q_object);
            % dq''/dphiop = d(R_str*q')/dthetaop  
            dqpp_dphiop = bsxpagemult(Rot_str_diff_phi, unit_q_object);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % q'' - Coordinates of each measurement in the SH centered coordinates
            x_sh_cut = q_pp(1, :, :);
            y_sh_cut = q_pp(2, :, :);
            q_pp = []; % free up memory
            % sin(theta'')
            sin_theta_sh_cut = sqrt(1-cos_theta_sh_cut.^2);
            % phi''
            phi_sh_cut = atan2(y_sh_cut, x_sh_cut);
            
            % [1 numofseg numOfvoxels]
            const_val = exp(1i.*phi_sh_cut) .* sin_theta_sh_cut .* cos_theta_sh_cut;
            
            % [numorders numofseg numOfvoxels]
            % block2 = 1   for order = 0
            %        = exp(i*phi'')*sin(theta'')*cos(theta'')*cos^[2*order](theta'')    for order > 2, 4, 6    .
            block2 = double([ones_struct; repmat(const_val, numOfCoeffs-1, 1, 1) .* block_cos_theta_powers(1:numOfCoeffs-1, :, :)]);
            const_val = [];
            
            % [] = sum_lm((  { alm*sqrt[(l-m)*(l+m+1)] * Ylm+1coeff } * exp(i*phi'')*sin(theta'')*cos(theta'')*cos^[2*order](theta'') ))
            % [1 numsegs numofvoxels] = [1 numorders numvoxels]*[numorders numofseg numOfvoxels]
            YlmPLUS1 = bsxpagemult(double(a_temp2), block2);
            block2 = []; % free up memory
            
            if any(m) % full calculations
                error('Calculations for m ~= 0 need to be checked')
                Ylm_temp = bsxpagemult(a_temp1, block_cos_theta_powers);
                
                %auxillary variables for gradient of error matrix in respect to
                %theta_struct and phi_struct in each voxel
                temp_theta_1 = dqpp_dthetaop(2, :, :) .* x_sh_cut - dqpp_dthetaop(1, :, :) .* y_sh_cut ./ x_sh_cut.^2;
                temp_theta_2 = dqpp_dthetaop(3, :, :);
                
                temp_phi_1 = dqpp_dphiop(2, :, :) .* x_sh_cut - dqpp_dphiop(1, :, :) .* y_sh_cut ./ x_sh_cut.^2;
                temp_phi_2 = dqpp_dphiop(3, :, :);
                
                expr1 = Ylm_temp .* (cos(phi_sh_cut).^2)*1i;
                expr2 = Ylm_temp .* (cos_theta_sh_cut./sin_theta_sh_cut.^2) + ...
                    YlmPLUS1 .* (exp(-1i*phi_sh_cut)./sin_theta_sh_cut);
                
                YlmPLUS1 = []; % free up memory
                phi_sh_cut = [];
                sin_theta_sh_cut = [];
                block_cos_theta_powers = [];
                
                %auxillary variable for the gradient (numOfvoxels x numofsegments)
                grad_theta_struct_aux_vol = real((expr1 .* temp_theta_1 - expr2 .* temp_theta_2) .* conj(sumlm_alm_Ylm));
                
                %auxillary variable for the gradient (numOfvoxels x numofsegments)
                grad_phi_struct_aux_vol = real((expr1 .* temp_phi_1 - expr2 .* temp_phi_2) .* conj(sumlm_alm_Ylm));
                
                expr1 = [];
                expr2 = [];
                
            else % optimized version - only valid for m = 0
                %auxillary variables for gradient of error matrix with respect to
                %theta_struct and phi_struct in each voxel
                % - [[[ sum_lm((  { alm*sqrt[(l-m)*(l+m+1)] * Ylm+1 } * exp(i*phi'')*sin(theta'')*cos(theta'')*cos^[2*order](theta'') )) ]]] * exp(-1i*phi'') ./ sin_theta'' .* conj( sum_lm(a_lm*Y_lm) )               
                % simplifying a bit     - [[[ sum_lm((  { alm*sqrt[(l-m)*(l+m+1)] * Ylm+1 } * cos(theta'')*cos^[2*order](theta'') )) ]]] * conj( sum_lm(a_lm*Y_lm) )  
                % MGS - note that here for general m is already missing the
                % first part of Eq(12a) of Liebi 2018
                expr2 = - YlmPLUS1 .* exp(-1i*phi_sh_cut) ./ sin_theta_sh_cut .* conj(sumlm_alm_Ylm);
                
                YlmPLUS1 = []; % free up memory
                phi_sh_cut = [];
                sin_theta_sh_cut = [];
                block_cos_theta_powers = [];
                
                %auxillary variable for the gradient (numOfvoxels x numofsegments)
                %   = expr2 * dq''/dthetaop 
                %   = - dq''/dthetaop * [[[ sum_lm((  { alm*sqrt[(l-m)*(l+m+1)] * Ylm+1coeff } * cos(theta'')*cos^[2*order](theta'') )) ]]] * conj( sum_lm(a_lm*Y_lm) )                                                 .
                %   = - dq''/dthetaop * [[[ sum_lm((   alm*  { sqrt[(l-m)*(l+m+1)] * {{Ylm+1coeff  * cos(theta'')*cos^[2*order](theta'')}} )) ]]] * conj( sum_lm(a_lm*Y_lm) )                                                 .
                % MGS - note, this expression is later multiplied by the
                % backprojected poisson difference (aux_grad_poisson_vol),
                % which is the first two lines of Eq(11) in Liebi et al
                % 2018. 
                % Line 3 of that equation is [dq''/dthetaop]*[dYlm/dtheta],
                % the second part is Eq. (12a). When comparing the first
                % part to the calculation above I see two things, one is
                % the missing exp(-iphi) but this dependance cancels from
                % Ylm+1 function, I also see one too many alm and one too many cos(theta''). I also dont
                % know why I have the first sum over lm, this , maybe this
                % sum is missing in the paper? or maybe I messed up above.
                % Line 4 is equal to zero because m = 0, Eq. (12b)
                % Then the last line of that equation is conj(
                % sum_lm(a_lm*Y_lm) ).
                grad_theta_struct_aux_vol = real(expr2 .* dqpp_dthetaop(3, :, :));
                
                %auxillary variable for the gradient (numOfvoxels x numofsegments)
                %   = - dq''/dphiop * [[[ sum_lm((  { alm*sqrt[(l-m)*(l+m+1)] * Ylm+1 } * cos(theta'')*cos^[2*order](theta'') )) ]]] * conj( sum_lm(a_lm*Y_lm) )                                                 .
                grad_phi_struct_aux_vol = real(expr2 .* dqpp_dphiop(3, :, :));
                
                expr2 = [];
            end
            
            dqpp_dthetaop = [];
            dqpp_dphiop = [];
            
            grad_theta_struct_aux_vol = permute(grad_theta_struct_aux_vol, [3, 2, 1]);
            grad_phi_struct_aux_vol = permute(grad_phi_struct_aux_vol, [3, 2, 1]);
            
            grad_theta_struct_aux_vol = reshape(grad_theta_struct_aux_vol, ny, nx, nz, numOfsegments);
            grad_phi_struct_aux_vol   = reshape(grad_phi_struct_aux_vol, ny, nx, nz, numOfsegments);
            
            %%%%%gradient of error matrix in respect to theta_struct and phi_struct, with
            %%%  error metric e=2*(sqrt(I(estimate))-sqrt(I(data)))^2) poisson noise
            %%%  and data_synt=abs((Ymn_2D).^2), for general m
            grad_theta_struct = grad_theta_struct + 4*sum(aux_grad_poisson_vol.*grad_theta_struct_aux_vol, 4);
            grad_phi_struct = grad_phi_struct + 4*sum(aux_grad_poisson_vol.*grad_phi_struct_aux_vol, 4);
            
            aux_grad_poisson_vol = []; % free up memory
            grad_theta_struct_aux_vol = [];
            grad_phi_struct_aux_vol = [];
        end
    end
    
    sumlm_alm_Ylm = []; % free up memory
end

E_temp = E;
if find_grad
    e = optimization.errorplot(E);
    iteration = length(e);
    fprintf('*************************************************** \nIteration %d \n', iteration);
else
    fprintf('Line Search \n')
end

if return_Ereg || p.regularization_angle
    if find_grad
        fprintf('The regularization parameter is %.3d \n', p.regularization_angle_coeff)
        [Reg, grad_theta_reg, grad_phi_reg] = angle_regularization(...
            theta_struct, phi_struct, sin_theta_struct, cos_theta_struct, sin_phi_struct, cos_phi_struct, ...
            ny, nx, nz);
        grad_theta_struct = grad_theta_struct + p.regularization_angle_coeff.*grad_theta_reg./numOfvoxels;
        grad_phi_struct = grad_phi_struct + p.regularization_angle_coeff.*grad_phi_reg./numOfvoxels;
    else
        Reg = angle_regularization(...
            theta_struct, phi_struct, sin_theta_struct, cos_theta_struct, sin_phi_struct, cos_phi_struct, ...
            ny, nx, nz);
    end
    
    Ereg = (p.regularization_angle_coeff.*Reg)/numOfvoxels;
    
    E = E + Ereg;
    fprintf('          total_error = data_error + regul_error \n');
    fprintf('          %f   = %f  + %f \n', E, E_temp, Ereg);
else
    fprintf('          data_error = %f \n', E_temp);
end

if find_grad
    %sieves regularization, blurring the gradient
    if find_coefficients && p.regularization
        if p.slice
            grad_a(:,:,p.slice_pos,:) = convn(grad_a(:,:,p.slice_pos,:), p.kernel, 'same');
        else
            grad_a = convn(grad_a, p.kernel, 'same');
        end
    end
    
    % mask gradients
    if find_coefficients
        grad_a = grad_a .* mask3D;
    end
    
    if find_orientation
        grad_theta_struct = grad_theta_struct .* mask3D;
        grad_phi_struct = grad_phi_struct .* mask3D;
    end
end

if find_coefficients
    a_tomo = reshape(a, ny, nx, nz, numOfCoeffs);
end

if find_grad
    % plotting
    if (p.display ~= 0) && (mod(iteration, p.dispinter) == 0)
        if p.plot
            figure(102)
            hold on
            subplot(2, 2, 1)
            
            if find_coefficients
                imagesc(sum(grad_a(:,:,:,1),3)); %plot the first coefficent
                title(sprintf('gradient of a0,2 \n data_error = %0.2d ', E_temp),  'interpreter','none')
            else
                imagesc(real(sum(grad_theta_struct,3))); %plot the first coefficent
                title(sprintf('Polar angle gradient \n data_error = %0.3f ',  E_temp),  'interpreter','none')
            end
            
            axis xy equal tight
            colorbar
            colormap jet
            drawnow
            subplot(2, 2, 2)
            
            if find_coefficients
                if numOfCoeffs > 1
                    t =(sum(((a_tomo.*optimization.spherical_harm(0,0,0,0)).^2),4)); %plot the first coefficent
                    imagesc(sum(t, 3))
                else
                    imagesc(sum(((a_tomo.*optimization.spherical_harm(0,0,0,0)).^2),3)); %plot the first coefficent
                end
                title(sprintf('Calculated SH \n data_error = %0.2d ',   E_temp),  'interpreter','none')
            else
                imagesc(real(sum(grad_phi_struct,3))); %plot the first coefficent
                title(sprintf('Azimuthal angle gradient \n data_error = %0.3f ',   E_temp), 'interpreter','none')
            end
            
            axis xy equal tight
            colorbar
            colormap jet
            drawnow
            subplot(2, 2,[3, 4])
            plot(iteration, E_temp, 'o--b');
            
            if find_orientation
                if p.regularization_angle
                    hold on
                    loglog(iteration, Ereg, 'o--r');
                    hold off
                    legend({'data_error', 'regul_error'},'interpreter','none');
                else
                    legend({'data_error'},'interpreter','none');
                end
                
                if find_coefficients
                    title('Angles and coefficients optimization')
                else
                    title('Angles optimization')
                end
            else
                title('Coefficient optimization')
            end
            
            xlabel('Iteration number');
            ylabel('Error');
            grid on
            
        end
    end
end

% Apply soft limits to the coefficients 
if find_coefficients
    for ii = 1:numOfCoeffs
        a_slice = a_tomo(:,:,:,ii);

        ind_low = a_slice < p.coeff_soft_limits_low(ii);
        ind_high = a_slice > p.coeff_soft_limits_high(ii);

        a_slice(ind_low) = a_slice(ind_low) - p.coeff_soft_limits_low(ii);
        a_slice(ind_high) = a_slice(ind_high) - p.coeff_soft_limits_high(ii);
        a_slice(~ind_low & ~ind_high) = 0;

        E = E + p.soft_limit_weight_coeff(ii)*sum(a_slice(:).^2);

        if find_grad
            grad_a(:,:,:,ii) = grad_a(:,:,:,ii) + 2*p.soft_limit_weight_coeff(ii)*a_slice;
        end
    end
end

if find_grad
    % gradient output vector
    grad = [];
    if find_orientation
        grad_theta_struct = reshape(grad_theta_struct, 1, []);
        grad_phi_struct = reshape(grad_phi_struct, 1, []);
        grad = [grad, grad_theta_struct, grad_phi_struct];
    end

    % Remove from tne gradient the coefficients that are not optimized
    % MGS Notice that here we miss an opportunity for optimization,
    % all the gradients are calculated above, even those for orders that
    % will not be used.
    %RSD: Do not understand this part entirely... Is code not optimised?
    %RSD: Believe is unnecessary. 
    if find_coefficients
        grad_a_temp = zeros(ny, nx, nz, numOfordersopt);
        counter = 1;
        for ii = 1:numel(p.opt_coeff)
            if p.opt_coeff(ii)
                grad_a_temp(:,:,:,counter) = grad_a(:,:,:,ii);
                counter = counter + 1;
            end
        end
        grad_a = reshape(grad_a_temp, 1, []);
        grad = [grad, grad_a];
    end

    % saving intermediate results
    filepath = p.optimization_output_path;
    if find_orientation && find_coefficients
        save([filepath 'current_opt_saxs_all.mat'], 'a_tomo', 'grad_a', ...
            'grad_theta_struct','grad_phi_struct','theta_struct','phi_struct','e');
        
    elseif find_orientation && ~find_coefficients
        save([filepath 'current_opt_saxs_angle.mat'], ...
            'grad_theta_struct','grad_phi_struct','theta_struct','phi_struct','e');
        
    elseif ~find_orientation && find_coefficients
        if numOfCoeffs == 1
            save([filepath 'current_opt_symsaxs.mat'], 'a_tomo', 'grad_a', 'e');
        else
            save([filepath 'current_opt_saxs.mat'], 'a_tomo', 'grad_a', 'e');
        end
    end
end



