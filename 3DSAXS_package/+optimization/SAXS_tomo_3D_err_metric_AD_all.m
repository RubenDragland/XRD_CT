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

%RSD: Memory management.
%aux_vars = {'x','y','z' };
%clear(aux_vars{:});


%%% RSD: Change of order finished

%RSD: Implement if elseif else with the options being optimise all, orientation, coefficients, or not find grad
if p.python && find_grad && p.batch
    
    % RSD: Make path automatic. 
    P = py.sys.path;
    search_path = [p.parent '/Autodiff_package'];
    if count(P,search_path) == 0
        insert(P,int32(0),search_path);
    end
    
    py.importlib.import_module("batch_AD");
    
    theta_struct_it = py.numpy.array(theta_struct);
    phi_struct_it = py.numpy.array(phi_struct);
    a_temp_it = py.numpy.array(a_temp);
    Xi = py.numpy.array(X);
    Yi = py.numpy.array(Y);
    Zi = py.numpy.array(Z);
    unit_q_beamline_i = py.numpy.array(unit_q_beamline);
    Ylm_coef_i = py.numpy.array(Ylm_coef);
    
    fields = {'diode', 'par', 'fnames', 'integ'};
    %batch_projection = rmfield(projection, fields);
    batch_projection = 1; %p.projection_filename; % batch_projection(1);
    
    batch_results = py.batch_AD.SH_main(theta_struct_it, phi_struct_it, a_temp_it, ny, nx, nz, numOfsegments, batch_projection, p, Xi, Yi, Zi, numOfpixels, unit_q_beamline_i, Ylm_coef_i, find_coefficients, find_orientation, numOfCoeffs, numOfvoxels, find_grad);
    
    E = batch_results{1};
    
    if ~isempty( double(batch_results{2}) )
        AD_grad_coeff = double(batch_results{2});
        grad_a = grad_a + reshape( permute(AD_grad_coeff, [3,2,1]) , ny, nx, nz, numOfCoeffs) ; 
     end
    if ~isempty( double(batch_results{3}) )
        AD_grad_theta = double( batch_results{3});
        AD_grad_phi = double (batch_results{4} );
        grad_theta_struct = grad_theta_struct + reshape( permute(AD_grad_theta, [3,2,1]), ny, nx, nz); 
        grad_phi_struct = grad_phi_struct + reshape( permute(AD_grad_phi, [3,2,1]), ny,nx,nz);
    end
    
    
elseif p.python && find_grad
    
    P = py.sys.path;
    if count(P,'C:\Users\Bruker\OneDrive\Dokumenter\NTNU\XRD_CT\Autodiff_package') == 0
        insert(P,int32(0),'C:\Users\Bruker\OneDrive\Dokumenter\NTNU\XRD_CT\Autodiff_package');
    end
    
    py.importlib.import_module("forward_backward_AD");
    
    for ii = 1:length(projection)
    %parfor ii = 1:length(projection) % RSD: Set variables before used.
    %Cannot use parfoor. Make the python version entirely batched instead. 

        theta_struct_it = py.numpy.array(theta_struct);
        phi_struct_it = py.numpy.array(phi_struct);
        a_temp_it = py.numpy.array(a_temp);
        Xi = py.numpy.array(X);
        Yi = py.numpy.array(Y);
        Zi = py.numpy.array(Z);
        unit_q_beamline_i = py.numpy.array(unit_q_beamline);
        Ylm_coef_i = py.numpy.array(Ylm_coef);

        current_projection = projection(ii);

        %[error_norm, AD_grad_coeff, AD_grad_theta, AD_grad_phi] = py.forward_backward_AD.main(pyargs("theta_struct_it",theta_struct_it, "phi_struct_it", phi_struct_it, "a_temp_it", a_temp_it,"ny", ny,"nx", nx, "nz", nz, "numOfsegments", numOfsegments,"current_projection", current_projection, "p", p, "X", X, "Y", Y,"Z", Z, "numOfpixels", numOfpixels, "unit_q_beamline", unit_q_beamline, "Ylm_coef", Ylm_coef, "find_coefficients", find_coefficients, "find_orientation", find_orientation, "numOfCoeffs", numOfCoeffs, "numOfvoxels", numOfvoxels) ); %py.forward_backward_AD.main(theta_struct_it, phi_struct_it, a_temp_it, ny, nx, nz, numOfsegments, current_projection, p, X, Y, Z, numOfpixels, unit_q_beamline, Ylm_coef, find_coefficients, find_orientation, numOfCoeffs, numOfvoxels);
        %res = py.forward_backward_AD.print_test("abc");
        %[error_norm, AD_grad_coeff, AD_grad_theta, AD_grad_phi] = py.forward_backward_AD.main(theta_struct_it, phi_struct_it, a_temp_it, ny, nx, nz, numOfsegments, current_projection, p, X, Y, Z, numOfpixels, unit_q_beamline, Ylm_coef, find_coefficients, find_orientation, numOfCoeffs, numOfvoxels);
        imeres =  py.forward_backward_AD.main(theta_struct_it, phi_struct_it, a_temp_it, ny, nx, nz, numOfsegments, current_projection, p, Xi, Yi, Zi, numOfpixels, unit_q_beamline_i, Ylm_coef_i, find_coefficients, find_orientation, numOfCoeffs, numOfvoxels, find_grad) ;
        
        %[error_norm, AD_grad_coeff, AD_grad_theta, AD_grad_phi] = [imeres{1}, imeres{2}, imeres{3}, imeres{4}];
        
        error_norm = imeres{1};
        E = E + error_norm;
        if ~isempty( double(imeres{2}) )
            AD_grad_coeff = double(imeres{2});
            grad_a = grad_a + reshape( permute(AD_grad_coeff, [3,2,1]) , ny, nx, nz, numOfCoeffs) ; % RSD: Check this!
            %RSD: Hacking alarm: This gives the correct output format given
            %RSD: the mix of C- and F- reshaping inside the python code.
            %RSD: Status quo everything is still shady. 
        end
        if ~isempty( double(imeres{3}) )
            AD_grad_theta = double( imeres{3});
            AD_grad_phi = double (imeres{4} );
            grad_theta_struct = grad_theta_struct + reshape( permute(AD_grad_theta, [3,2,1]), ny, nx, nz); % RSD: Probably fix this
            grad_phi_struct = grad_phi_struct + reshape( permute(AD_grad_phi, [3,2,1]), ny,nx,nz);
        end

    end

elseif find_grad && find_orientation && find_coefficients

    [E, grad_a, grad_theta_struct, grad_phi_struct] = SAXS_AD_all_forward_backward(theta_struct, phi_struct, a_temp, ny, nx, nz, numOfsegments, projection, p, X, Y, Z, numOfpixels, unit_q_beamline, Ylm_coef, find_coefficients, numOfCoeffs, numOfvoxels); %RSD: Placeholder

elseif find_grad && find_orientation

    [E, ~, grad_theta_struct, grad_phi_struct] = SAXS_AD_all_forward_backward(theta_struct, phi_struct, a_temp, ny, nx, nz, numOfsegments, projection, p, X, Y, Z, numOfpixels, unit_q_beamline, Ylm_coef, find_coefficients, numOfCoeffs, numOfvoxels); %RSD: Placeholder. Idea, only one file with some if/else. 

%elseif find_grad && find_coefficients This part we already have. grad and no grad is finished implemented within each other.

else

    %RSD: This code calculates the coefficient gradients or no gradients. Move most of the general code outside of the if-else.

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

    %RSD: If no regularization, consider to delete sin_theta_struct etc after creating the rot_str matrix to save memory. 
    %aux_vars = {"sin_theta_struct", "cos_theta_struct", "sin_phi_struct", "cos_phi_struct" };
    %clear(aux_vars{:})

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
        %   [numOforder numOfsegments numOfvoxels] = [numOforder numOforder]*[numOforder numOfsegments numOfvoxels]    .
        Ylm = bsxpagemult(double(Ylm_coef), block_cos_theta_powers);  % RSD: RETRIEVE THE SH FUNCTION FROM MATRIX PRODUCT WITH ARRAY OF POWERS OF COSINE. YLM-COEFFICIENTS MATCHED WITH ITS "COSINE-POLYNOMIAL"

        %aux_vars = {"q_pp", "cos_theta_sh_cut", "unit_q_object", "block_theta_powers"};
        %clear(aux_vars{:})
        
        % sum_lm(a_lm*Y_lm)
        %   [1 numOfsegments numOfvoxels] = [1 numOforders numOfvoxels] * [numOforder numOfsegments numOfvoxels]
        
        cost_function = @SAXS_AD_cost_function;
        current_projection = projection(ii) ; % RSD: Only broadcast the current projection to avoid parfor crash. 
        [error_norm, AD_grad_coeff] = SAXS_AD_forward_backward(cost_function, a_temp, Ylm, ny, nx, nz, numOfsegments, data, current_projection, Rot_exp_now, p, find_grad, X, Y, Z, numOfpixels);
        
        if find_grad
            error_norm = extractdata(error_norm); 
        end
        E = E + error_norm;

        %%%the gradients RSD: Many changes to the code here to finish the automatic differentiation.
        if find_grad
            
            if find_coefficients % RSD: May be abundant if, but keep for now
                
                %RSD: Need to reshape
                AD_grad_coeff = extractdata(AD_grad_coeff); % RSD: No longer need for dlarrays
                grad_a = grad_a + reshape( permute(AD_grad_coeff, [3,2,1]), ny,nx,nz,numOfCoeffs); %permute( reshape(AD_grad_coeff, ny, nx, nz, numOfCoeffs), [ %reshape(AD_grad_coeff, ny, nx, nz, numOfCoeffs); %RSD: MATLAB reshape fucking shit...
            end %RSD: Since a_temp has shape 1xnumCoeffsx(nynznz), permutation before reshape is necessary. 
            Ylm = []; %#ok<*NASGU> % free up memory
        end
    end

end            %RSD: Or possibly escape the if-statement here.
            



E_temp = E;
if find_grad
    global Err_hist;
    Err_hist = [Err_hist, E];
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



