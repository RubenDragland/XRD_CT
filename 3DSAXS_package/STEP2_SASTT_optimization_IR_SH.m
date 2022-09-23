% WIP: Optimization with the Ir method as initial guess

%*-------------------------------------------------------------------------------------*
%|                                                                                     |
%|  Except where otherwise noted, this work is licensed under a                        |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0                          |
%|  International (CC BY-NC-SA 4.0) license.                                           |
%|                                                                                     |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)                  |
%|                                                                                     |
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

error('Please run this script by section');

%% Load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent = cd ;
base_path = '/Data sets/' ; %'~/Data10'; % = ~/Data10 for online analysis, provide the path for offline analysis Ex: '/das/work/p16/p16649/'
base_path = [parent base_path] ; % SAFE WITH FULL PATH
sample_name = 'SASTT_carbon_knot_aligned_ASTRA_corrected' ; %'sample_name'; % name given in the saxs_caller_template
add_name = '';       % additional name the optimizations: = [ ] if not needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load data
%filename = fullfile(base_path,sprintf('analysis/SASTT/%s/projection_data/SASTT_%s_aligned_ASTRA.mat', ...sample_name, sample_name));
filename = sprintf('%s%s.mat', base_path, sample_name); %RSD: CURRENTLY WHAT IS NEEDED.
load(filename);

params_IRTT.results = fullfile(base_path, sprintf('analysis/SASTT/%s/IRTT/%s/results/', sample_name, add_name)); % MAY ACTUALLY BE A GOOD IDEA... KEEP /OUTPUT/ EMPTY FOR NOW
if ~exist(params_IRTT.results, 'file')
    mkdir(params_IRTT.results);
end
params_IRTT.figures = fullfile(base_path, sprintf('analysis/SASTT/%s/IRTT/%s/figures/', sample_name, add_name));
if ~exist(params_IRTT.figures, 'file')
    mkdir(params_IRTT.figures);
end

%% Visualize aligned segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stepping = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure (12)
for ii = 1:stepping:length(projection)
    imagesc(mean(projection(ii).data, 3).*projection(ii).window_mask)
    title(sprintf('Projection %d/%d', ii, length(projection)))
    axis equal tight xy
    drawnow
end

%% Tensor reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct tensor model from projection of all segments.
% Set a lower iteration first to check the result, if it's correct then
% run more iterations to minimize the error.
num_iter = 10000;
% Typical reconstruction would need 1000 iterations to see the rough
% model, and 10000 iterations for final result.
p.skip_projections = 1;
if_show = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parallel_scattering = 1; % Parallel scattering must be set to 1 if used for initial guess of SH reconstruction.

[tomotensor, B_segs] = tensor_reconstruct(projection(1:p.skip_projections:length(projection)), ...
    num_iter, if_show,[],[],parallel_scattering);

% save figure
filename = sprintf('%sIRTT_%s_q_%d-%d_%s', params_IRTT.figures, sample_name, projection(1).par.r_sum{1}(1),projection(1).par.r_sum{1}(end), add_name);
        print(gcf, filename, '-dpng','-r300');

% The format of tomotensor output is a  [Z, X, Y, T] array of a symmetric
% tensor for each voxel in the model.
% For the last dimension T the elements are arranged as (Txx,Tyy,Tzz,Txy,Txz,Tyz)

%% SAVE the tensor model output (eigen vector map) and export to paraview
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visulizing the tensor directly is not easy, therefore here it solves a
% certain eigen vector of the tensor and visualize the vector map.
threshold_low = 0;                          % The threshold for the eigenvalue to remove background.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the three eigen vectors
tomotensor_vecs_1 = solve_eig(tomotensor, 1, threshold_low);
tomotensor_vecs_2 = solve_eig(tomotensor, 2, threshold_low);
tomotensor_vecs_3 = solve_eig(tomotensor, 3, threshold_low);

% save the results
% next tensor tomo
filename = sprintf('%sIRTT_%s_q_%d-%d_%s.mat', params_IRTT.results, sample_name, projection(1).par.r_sum{1}(1),projection(1).par.r_sum{1}(end), add_name);
save(filename, 'tomotensor', 'tomotensor_vecs_1', 'tomotensor_vecs_2', ...
    'tomotensor_vecs_3', 'B_segs', 'num_iter', 'parallel_scattering', '-v6');

% To Export: paraview (ellipsoids)
fid = fopen(sprintf('%s/IRTT_paraview_data_%s_%s.raw', params_IRTT.results, ...
    sample_name,add_name), 'w');
data_paraview = permute(tomotensor_vecs_1, [4, 2, 1, 3]);
fwrite(fid,data_paraview,'float');
fclose(fid);

%% Create 3D mask from IRTT tomogram
% input is for the spherical harmonics optimization is  Y X Z, but tomotensor
% is [X, Y, Z, 6]
tomotensor_vecs_1_abs = sum(tomotensor_vecs_1.^2,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
make_3Dmask = 1; % 1 to create a 3D mask, 0 to see the symmetric intensity tomogram
mask.cut_off = 1; % = [] for automatic detection of the mean value using Otsu threshold
mask.diam_cylinder = 50; % Diameter of a cylinder outside, = [] does not create the mask
mask.gauss_filter = 1; % Apply a filter to the data to reduce noise
mask.openning = 1;     % Apply an kind of openning operation on the mask, this makes it larger and is important to not have too tight constraint.
coloraxis = [0, 500]; % apply same color scale to the three projections xy, yz and xz fron view3
viewpar.interactive = 1;  % = true if you want to explore the volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if make_3Dmask
   [s.mask3D, cut_off] = create_mask3D(tomotensor_vecs_1_abs, mask);
else
    [s.mask3D, ~] = ones(size(tomotensor_vecs_1_abs));
    cut_off = 0;
end

% plot the masked symmetriuc intensity reconstruction
view3((s.mask3D.*tomotensor_vecs_1_abs - cut_off*3*(1-s.mask3D)), coloraxis,viewpar)

%% Give as input for SH optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimize all coefficients (a0, a2, a4 and a6) and angles (theta and phi)
% whihout any constriction
apply_angle_regularization_SH = 0; %RSD: 1 ERROR BECAUSE ARRAY DOES NOT EXIST YET... ; % =1 to regularize angles, =0 not to
p.skip_projections = 1;              % to reduce the number of projections

p.opt_coeff = [1, 1, 1, 1];

p.find_orientation = 1;      % Optimize over the main orientation of the structure
p.find_coefficients = any(p.opt_coeff);     % Optimize over coefficients
p.regularization = 0;        % Sieves regularization on coeff

%parameters for the sieves regularization (blurring the gradient)
kernel3D = window3(5,5,5,@hamming);
p.kernel = kernel3D./sum(kernel3D(:)); % for normalization (sum equals 1)

p.itmax = 50;   % maximum number of iterations: about 50-100

p.avoid_wrapping = 1;     % avoid wrapping over 2Pi of the angle

% HAS TO BE A SMALL ERROR IN THE CODE. P.RESULTS DOES NOT EXIST. COPY IT FROM IRTT.RESULTS
p.results = params_IRTT.results ;
p.sample_name = sample_name ;
p.add_name = add_name ;
p.figures = params_IRTT.figures ;
all_again_filename = sprintf('%s/result_%s_q%d-%d_all_again_%s.mat', p.results, ...
    p.sample_name, projection(1).par.r_sum{1}(1),  projection(1).par.r_sum{1}(end), p.add_name);

p.save.output_filename = all_again_filename;
p.save.image_filename = sprintf('%s/optimization_all_q%d-%d_%s', p.figures, ...
    projection(1).par.r_sum{1}(1),  projection(1).par.r_sum{1}(end), p.add_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the IRTT results
initial_guess = sprintf('%sIRTT_%s_q_%d-%d_%s.mat', params_IRTT.results, sample_name, projection(1).par.r_sum{1}(1),projection(1).par.r_sum{1}(end), add_name);
IG = load(initial_guess);

if ~isfield(IG,'parallel_scattering')
    fprintf('WARNING: cannot determine scattering direction used in tensor reconstruction.\n');
elseif IG.parallel_scattering==0
    fprintf('WARNING: WRONG scattering direction used in tensor reconstruction.\n');
end

[IG.theta_out1, IG.phi_out1, IG.a_out1] = Tensor_to_SH(IG.tomotensor);

p.points_per_line = size(tomotensor_vecs_1_abs, 1);
p.lines = size(tomotensor_vecs_1_abs, 2);
p.numOfvoxels = p.points_per_line*p.lines*p.lines;

% prepare params for SH
s.theta.data = IG.theta_out1;
s.phi.data = IG.phi_out1;
for ii = 1:numel(p.opt_coeff)
    if p.opt_coeff(ii)
        s.a(ii).data = IG.a_out1(1:p.numOfvoxels);
        IG.a_out1 = IG.a_out1(p.numOfvoxels+1:end);
    end
end

% RSD: EDIT: SOME DIFFERENCES BETWEEN IR_SH PREPARATIONS AND SH_PREPARATIONS. RESULTS IN ERRORS. HENCE, INCLUDE SH_PREPARATIONS:
E = []; % X-RAY SOURCE ENERGY. THINK SHOULD BE EMPTY FOR SAXS
if isempty(E)
    projection(1).integ.theta_det = pi/2;
end
% NOTE THAT ONLY EMPTY E WILL NOW NOT RESULT IN ERROR. DO NOT UNDERSTAND THIS PART.

% END OF CODE EDIT
p.add_name = add_name;
p.sample_name = sample_name;

% prepare paths for saving the data
p.optimization_output_path = fullfile(base_path,sprintf('analysis/SASTT/%s/SH/%s/optimization_output/', sample_name, p.add_name));
if ~exist(p.optimization_output_path, 'file')
    mkdir(p.optimization_output_path);
end
p.figures = fullfile(base_path,sprintf('analysis/SASTT/%s/SH/%s/figures/', sample_name, p.add_name));
if ~exist(p.figures, 'file')
    mkdir(p.figures);
end
p.results = fullfile(base_path,sprintf('analysis/SASTT/%s/SH/%s/results/', sample_name, p.add_name));
if ~exist(p.results, 'file')
    mkdir(p.results);
end

% define regularization of angles
if apply_angle_regularization_SH
    p.regularization_angle = 1; %regularization of angle
    % calculate estimated regularization coefficient
    a_tomo = reshape(IG.a_out1(1:p.numOfvoxels), p.points_per_line, p.lines, p.lines);
    t = sum(((a_tomo.*optimization.spherical_harm(0,0,0,0)).^2),4); %plot the first coefficent
    tt = (sum(t, 3).*projection(1).window_mask); %axis tight equal xy
    mu_guess = mean(mean(tt));
    p.regularization_angle_coeff = mu_guess;  %find the appropriate mu with L-curve from step 2.6
else
    p.regularization_angle = 0; %regularization of angle
end

% run optimization
close all

% parameters for optimization of angles with the 4 coefficients
l = [0 2 4 6];  % Polar order
m = [0 0 0 0];  % Azimuthal order

for ii = 1:numel(l)
    s.a(ii).l = l(ii);
    s.a(ii).m = m(ii);
end

% needed when IRTT results are given as initial guess
if ~isfield(p, 'phi_det')
    % %%% 2D plot characteristics (the data) %%%
    % Half of the Angular segments in integration
    numsegments = length(projection(1).integ.phi_det)/2;
    p.numsegments = numsegments;
    % read the phi_det from integ
    % consider only the first half because of the symmetry assumption for the
    % segments
    p.phi_det = deg2rad(projection(1).integ.phi_det(1:numsegments)); %in radians
    %%%%%%%%%%%%%%%%%%%%
end

%optimize
[p, s] = optimization.optimize_SH(projection, p, s);
