% function [p, s] = optimize_SH(projection, p, s)
% 
% Inputs:
%** projection  Projection data structure
%** p           Optimization parameter structure
%** s           Tensor in spherical harmonics representation. 
%       s.a(ii).l   iith coefficient polar order
%       s.a(ii).m   iith coefficient azimuthal order
%       s.a(ii).data(ny,nx,nz)  iith coefficient volume data in cSAXS 
%                               beamline coordinates, ny is up, nz is the  
%                               beam propagation direction, and nx is to 
%                               the right looking into the beam
%       s.theta.data and s.phi.data     Angles theta_op and phi_op that 
%           parametrize the local preferential orientation, i.e. the local
%           orientation of the spherical harmonic azimuth, with respect to  
%           the object coordinates. See Eq.(2) in [1].
%       s.mask3D    Binary volume mask that defines the regions that can be
%                   optimized. A support constraint. Defined as M in [1].
% 
% returns:
% ++ E          Error metric
% ++ grad       Gradient of E with respect to optimization variables
% ++ proj_out   Synthetic projections from the 3D data, also includes the
%               2D error map and the data error per projection
% ++ Ereg       Orientation regularization error metric
%
%   [1] Liebi et al., "Small-angle X-ray scattering tensor tomography: 
%       model of the three-dimensional reciprocal-space map, reconstruction 
%       algorithm and angular sampling requirements," Acta Cryst. 12–24 (2018). 

%*-------------------------------------------------------------------------------------*
%|                                                                                     |
%|  Except where otherwise noted, this work is licensed under a                        |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0                          |
%|  International (CC BY-NC-SA 4.0) license.                                           |
%|                                                                                     |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)                  |
%|                                                                                     |
%|      Author: CXS group, PSI                                                         |
%*------------------------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If this code, or subfunctions or parts of it, is used for research in a
%   publication or if it is fully or partially rewritten for another
%   computing language the authors and institution should be acknowledged
%   in written form and additionally you should cite:
%     M. Liebi, M. Georgiadis, A. Menzel, P. Schneider, J. Kohlbrecher,
%     O. Bunk, and M. Guizar-Sicairos, “Nanostructure surveys of
%     macroscopic specimens by small-angle scattering tensor tomography,”
%     Nature 527, 349-352 (2015).   (doi:10.1038/nature16056)
% and
%     M. Liebi, M. Georgiadis, J. Kohlbrecher, M. Holler, J. Raabe, I. 
%     Usov, A. Menzel, P. Schneider, O. Bunk and M. Guizar-Sicairos,
%     "Small-angle X-ray scattering tensor tomography: model of the 
%     three-dimensional reciprocal-space map, reconstruction algorithm 
%     and angular sampling requirements," Acta Cryst. A74, 12-24 (2018). 
%     (doi:10.1107/S205327331701614X)
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p, s] = optimize_SH(projection, p, s)
% Number of coefficients
order = numel([s.a.l]);

% set default parameters in case they are missing
if ~isfield(p,'avoid_wrapping')
    p.avoid_wrapping = 1;
    fprintf('p.avoid_wrapping is empty, setting to default value p.avoid_wrapping = %d \n', p.avoid_wrapping)
end
if ~isfield(p,'regularization')
    p.regularization = 0;
    fprintf('p.regularization is empty, setting to default value p.regularization = %d \n', p.regularization)
end
if ~isfield(s,'mask3D')
    s.mask3D = ones(p.points_per_line,p.lines,p.lines);
    fprintf('s.mask3D is empty, setting to all ones \n')
end
if ~isfield(p,'regularization_angle')
    p.regularization_angle = 0;
    fprintf('p.regularization_angle is empty, setting to default value p.regularization_angle = %d \n', p.regularization_angle)
end
if ~isfield(p,'coeff_soft_limits_low') || numel(p.coeff_soft_limits_low) ~= order
    p.coeff_soft_limits_low = -Inf * ones(order, 1);
    fprintf('p.coeff_soft_limits_low is empty, setting all to default value of -Inf \n');
else
    if numel(p.coeff_soft_limits_low) == 1
        p.coeff_soft_limits_low = p.coeff_soft_limits_low * ones(order, 1);
    end
end
if ~isfield(p,'coeff_soft_limits_high') || numel(p.coeff_soft_limits_high) ~= order
    p.coeff_soft_limits_high = Inf * ones(order, 1);
    fprintf('p.coeff_soft_limits_high is empty, setting all to default value of Inf \n')
else
    if numel(p.coeff_soft_limits_high) == 1
        p.coeff_soft_limits_high = p.coeff_soft_limits_high * ones(order, 1);
    end
end
if ~isfield(p,'soft_limit_weight_coeff') || numel(p.soft_limit_weight_coeff) ~= order
    p.soft_limit_weight_coeff = 1e5 * ones(order, 1);
    fprintf('p.soft_limit_weight_coeff is empty, setting to default value of 1e5 \n')
else
    if numel(p.soft_limit_weight_coeff) == 1
        p.soft_limit_weight_coeff = p.soft_limit_weight_coeff * ones(order, 1);
    end
end
if ~isfield(p,'display')
    p.display = 1;
    fprintf('p.display is empty, setting to default value p.display = %d \n', p.display)
end
if ~isfield(p,'dispinter')
    p.dispinter = 1;
    fprintf('p.dispinter is empty, setting to default value p.dispinter = %d \n', p.dispinter)
end
if ~isfield(p,'plot')
    p.plot = 1;
    fprintf('p.plot is empty, setting to default value p.plot = %d \n', p.plot)
end
if ~isfield(p,'volume_upsampling')
    p.volume_upsampling = 1;
    fprintf('p.volume_upsampling is empty, setting to default value p.volume_upsampling = %d \n', p.volume_upsampling)
end
if ~isfield(p,'method')
    p.method = 'bilinear';
    fprintf('p.method is empty, setting to default value p.method = %s \n', p.method)
end
if ~isfield(p,'filter_2D')
    p.filter_2D = 1;
    fprintf('p.filter_2D is empty, setting to default value p.filter_2D = %d \n', p.filter_2D)
end
if ~isfield(p,'slice')
    p.slice = 0;
    fprintf('p.slice is empty, setting to default value p.slice = %d \n', p.slice)
end
if ~isfield(p,'save') || ~isfield(p.save,'output_filename')
    p.save.output_filename = '';
    fprintf('p.save.output_filename is empty, not saving output files\n')
end
if ~isfield(p,'save') || ~isfield(p.save,'image_filename')
    p.save.image_filename = '';
    fprintf('p.save.image_filename is empty, not saving image files\n')
end

%%% Checks %%%
if numel(s.a)~= numel(p.opt_coeff)
   error('Number of coefficients, numel(s.a), and the number of optimization flags, numel(p.opt_coeff), do not match');
end

numOfvoxels = p.numOfvoxels;

% Arranging the input for optimization
opt_inputs = []; % Vector with parameters to optimize
if p.find_orientation
    opt_inputs = [s.theta.data(:); s.phi.data(:)];
end

for ii = 1:numel(p.opt_coeff)
    if p.opt_coeff(ii)
        opt_inputs = [opt_inputs; s.a(ii).data(:)];
    end
end
opt_inputs = opt_inputs.'; % RSD: Transpose https://se.mathworks.com/help/matlab/matlab_prog/matlab-operators-and-special-characters.html

Nvol = size(s.a(1).data);

opt_projection = projection(1:p.skip_projections:length(projection));

% optimization
tic
optimization.errorplot; %clears the persistent variable for error storage
%%%%%optimization with gradient
itmax = p.itmax; %maximum number of iteration (empty for default = 50)
ftol = [];  %relative function tolerance (empty for default = 1e-3)
xtol = [];  %absolute solution tolerance (empty for default = 1e-3)
opt_out = optimization.cgmin1('optimization.SAXS_tomo_3D_err_metric', opt_inputs, itmax, ftol, xtol, p, s, opt_projection);
timing = toc; % MUST BE WRONG IF TIMING...

% Rearranging the solution
if p.find_orientation
    s.theta.data = reshape(opt_out(1:numOfvoxels),Nvol);
    s.phi.data   = reshape(opt_out(numOfvoxels+1:2*numOfvoxels),Nvol);
    opt_out = opt_out(2*numOfvoxels+1:end);
end

for ii = 1:numel(p.opt_coeff)
    if p.opt_coeff(ii)
        s.a(ii).data = reshape(opt_out(1:numOfvoxels),Nvol);
        opt_out = opt_out(numOfvoxels+1:end);
    end
end

% Save relevant data in the output file
if p.save.output_filename
    e = optimization.errorplot([]);
    save(p.save.output_filename, 's', 'p', 'e', 'timing', '-v6'); % SOMETHING GOES WRONG HERE. TRY DIFFERENT VERSION, DIFFERENT LOCATION, DIFFERENT VARIABLES. SOMETHING WRONG WITH LOCATION. 
    %save(sprintf('C:\\Users\\Bruker\\OneDrive\\Dokumenter\\NTNU\\XRD_CT\\Data sets\\analysis\\SASTT\\SASTT_carbon_knot_aligned_ASTRA_corrected\\SH\\optimization_output\\result_%s%s_symSAXS.mat', p.sample_name, p.add_name ),'s', 'p', 'e', 'timing', '-v6' ); 
    %sym_opt_filename = sprintf(fullfile('%s','\result_%s%s_symSAXS.mat'), p.optimization_output_path, p.sample_name, p.add_name);
    %save(sym_opt_filename, 's', 'p', 'e', 'timing', '-v6');
end

% Save solution image
if p.save.image_filename
    print(gcf, p.save.image_filename, '-dpng', '-r300');
end
