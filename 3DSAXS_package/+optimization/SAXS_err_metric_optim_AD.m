function [E, grad, proj_out, Ereg] = SAXS_tomo_3D_err_metric_AD_all(opt_inputs, p, s, projection)
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

%RSD: The goal of this error metric is to reduce the impact of the dl-library by defining dlarrays once, and computing dlgradient once
%RSD: Hopefully, this is more efficient. 


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
%RSD: Run without parallel code for debugging.

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
    
    if find_grad 
        %p.method = 'nearest';       % RSD: Ensure nearest interpolation when calculating gradients. 
        p.volume_upsampling = 0;    % RSD: Same applies to upsampling
        %p.filter_2D = 0;            % RSD: Filter is not yet implemented for AD. RSD: EDIT: Ready for testing. 
    
    end

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
numOfsegments = p.numsegments;          % number of segements the saxs data is averaged (default is 16)
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
a_temp = reshape(a.', 1, numOfCoeffs, []);


%initalize arrays for the gradients
grad_a = zeros(ny, nx, nz, numOfCoeffs);

grad_theta_struct = zeros(ny, nx, nz);
grad_phi_struct   = zeros(ny, nx, nz);

%RSD: Up to this point kept constant just because. Could be improved, but change as little as possible.


%%% RSD: General features moved up in order to apply for all cases.

zeros_struct = zeros(1, 1, numOfvoxels);

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

ones_struct = ones(1, numOfsegments, numOfvoxels);

% optimization
E = 0;

% q (unit vector), 3D coordinates of the measurements in reciprocal space and in beamline coordinates. 
% Size 3x8 [xyz x numOfsegments]. The third dimension is actually only needed to account for Ewald sphere curvature
% correction which is included in theta_det
unit_q_beamline = [sin(theta_det).*cos(phi_det); sin(theta_det).*sin(phi_det);cos(theta_det)]; 

N = [ny nx nz]; % Size of tomogram  
x = (1:N(2)) - ceil(N(2)/2); %RSD: Is everything fine with this?
y = (1:N(1)) - ceil(N(1)/2);
z = (1:N(3)) - ceil(N(3)/2);
[X, Y, Z] = meshgrid(x, y, z); % (x,y,z) (2,1,3)
% If x is defined as the 3rd dimension N(3), then it should appear in the
% third position of the meshgrid command. If x is defined as the first
% dimension N(1), then it appears in the second position of the meshgrid,
% because meshgrid is (x,y,z) and size is (y,x,z)


%RSD: Divide into two parts in order to still be able to use bsxpagemult etc. 
if find_grad && find_orientation

    if find_coefficients
        %RSD: Introduce dlarrays for a_temp
        a_temp = dlarray(a_temp);
    end

    %RSD: Introduce dlarrays for theta_struct and phi_struct
    theta_struct = dlarray(theta_struct); %RSD: Unformatted
    phi_struct = dlarray(phi_struct);     %RSD: Unformatted

    [E, AD_grad_coeff, AD_grad_theta, AD_grad_phi] = dlfeval(@SAXS_AD_optim_fb_orientation, theta_struct, phi_struct, a_temp, ny, nx, nz, numOfsegments, projection, p, X, Y, Z, numOfpixels, unit_q_beamline, Ylm_coef, find_coefficients, numOfCoeffs, numOfvoxels );

    if find_coefficients
        grad_a = reshape( extractdata( AD_grad_coeff ), [ny, nx, nz, numOfCoeffs] );
    end

    grad_theta_struct = reshape(extractdata( AD_grad_theta ), [ny, nx, nz]);
    grad_phi_struct = reshape(extractdata(AD_grad_phi), [ny, nx, nz]);

    E = extractdata(E);

else
    %RSD: This code calculates the coefficient gradients or no gradients. Move most of the general code outside of the if-else.
    if find_grad && find_coefficients
        %RSD: Introduce dlarrays for a_temp
        a_temp = dlarray(a_temp);

        [E, AD_grad_coeff] = dlfeval(@SAXS_AD_optim_fb_coefficients, a_temp, theta_struct, phi_struct, ny, nx, nz, numOfsegments, projection, p, find_grad, X, Y, Z, numOfpixels, unit_q_beamline, Ylm_coef, find_coefficients, numOfCoeffs, numOfvoxels );

        grad_a = reshape( extractdata( AD_grad_coeff ), [ny, nx, nz, numOfCoeffs] );
        E = extractdata(E);
    else
        [E, AD_grad_coeff] = SAXS_AD_optim_fb_coefficients( a_temp, theta_struct, phi_struct, ny, nx, nz, numOfsegments, projection, p, find_grad, X, Y, Z, numOfpixels, unit_q_beamline, Ylm_coef, find_coefficients, numOfCoeffs, numOfvoxels );
    end
    
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



