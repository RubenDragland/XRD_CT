
%% Step 2.1a Load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent = cd;
base_path = '/Data sets/';  % = '~/Data10/' for online analysis,
base_path = [parent base_path] ;
                                           % provide the path for offline analysis
                                            % Ex: '/das/work/p16/p16649/'
sample_name = 'SASTT_carbon_knot_aligned_ASTRA_correctedny4nx4'; %'Validation_periodic_filter1_3cube_4off_0align_stripped'; % 'SASTT_carbon_knot_aligned_ASTRA_correctedny4nx4'; %'Synthetic_sample_ny4_ny4_all_coeffs'; %'sample_name';     % name given in the saxs_caller_template
p.add_name = '';        % additional name the optimizations: = [ ] if not needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../base
addpath utilities/
p.sample_name = sample_name;
% create folders to save data
% CHANGE HERE TO MOVE RESULTS TO DATA SETS FOLDER OUTSIDE OF PACKAGE
p.optimization_output_path = fullfile(base_path,sprintf('/analysis/SASTT/%s/SH/%s/optimization_output/', sample_name, p.add_name));
if ~exist(p.optimization_output_path, 'file')
    mkdir(p.optimization_output_path);
end
p.figures = fullfile(base_path,sprintf('/analysis/SASTT/%s/SH/%s/figures/', sample_name, p.add_name));
if ~exist(p.figures, 'file')
    mkdir(p.figures);
end
p.results = fullfile(base_path,sprintf('/analysis/SASTT/%s/SH/%s/results/', sample_name, p.add_name));
if ~exist(p.results, 'file')
    mkdir(p.results);
end


% load data
%filename = fullfile(base_path,sprintf('/analysis/SASTT/%s/projection_data/SASTT_%s_aligned_ASTRA.mat', ...
%    sample_name, sample_name));
filename = sprintf('%s%s.mat', base_path, sample_name);
fprintf('****** Step 2.1 Load the data ****** \n')
fprintf('Loading %s\n',filename)
load(filename);

p.projection_filename = filename;



%% Step 2.1b prepare SH structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E = [];%12.4;   % X-ray energy in keV, needed for Ewald sphere correction, leave 
            % empty and curvature of Ewald sphere will be ignored 
            % SHOULD BE EMPTY?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%% 2D plot characteristics (the data) %%%
% Half of the Angular segments in integration
p.numsegments = 8; % length(projection(1).integ.phi_det)/2; RSD: Another convention in the data sets provided by fredrik.

% read the phi_det from integ
% consider only the first half because of the symmetry assumption for the
% segments
p.phi_det = deg2rad(projection(1).integ.phi_det(1:p.numsegments)); %in radians

% for reconstruction of WAXS data (curvature of Ewald sphere)
% for SAXS put p.theta_det = [], this will apply the approximation of %DO NOT KNOW HOW TO HANDLE "E" AND THETA_DET...
% a planar cut (theta_det = pi/2);
if isempty(E)
    projection(1).integ.theta_det = pi/2;
    % RSD: SHOULD LINES 75-76 REGARDING R AND Q RANGE BE INCLUDED HERE?
else
    lambda = 1.23984197/E;  % wavelength of X-rays in nm
    %%% MGS - How general should this be? Should we allow each projection
    %%% to have different phi_det and theta_det?
    r_range = projection(1).par.r_sum{projection(1).par.which_qrange};
    q_range = projection(1).integ.q(r_range);
    %     q_value = projection(1).integ.q(1,projection(1).par.qresolved_q{1,params.qnow}(1));
    projection(1).integ.theta_det = pi/2 + mean( asin(q_range*lambda/(4*pi)) ); % for SAXS = pi/2, for WAXS = pi/2 - half the scattering angle;
end

% define the size of the tomogram (not the same as the tomogram, as some projection might change the size)
p.ny   = max(arrayfun(@(proj) size(proj.data, 1), projection));%+6
p.nx   = max(arrayfun(@(proj) size(proj.data, 2), projection));%+6
p.nz = p.nx;

% number of voxels in the tomogram
p.numOfvoxels = p.nx * p.ny * p.nz;

s = struct;

s.mask3D = ones(p.ny, p.nx, p.nz); %RSD: NOTE Y FIRST...

% parameters for optimization of 1 coefficient for symmetric intensity
% Spherical harmonic parameters for optimizing only the coefficient
% this is the same as fitting a voxel with a certain intensity - no angle
% available

l = 0;  % Polar order: degree
m = 0;  % Azimuthal order
a = 0.0001;   % initial guess for Coefficient

%%% Initial values for the SH angles
theta_init = pi/4; %45 degrees:  Polar orientation of structure
phi_init = pi/4; %45 degrees: Azimuthal orientation of structure

% calculate the initial values
s.theta.data = ones(p.ny,p.nx,p.nz)*theta_init;  % Polar orientation of structure
s.phi.data   = ones(p.ny,p.nx,p.nz)*phi_init;      % Azimuthal orientation of structure

global Err_hist;
Err_hist = [];

for ii = 1:numel(l)
    s.a(ii).data = ones(p.ny,p.nx,p.nz)*a(ii);
    s.a(ii).l = l(ii);
    s.a(ii).m = m(ii);
end


% Q-resolved %RSD: Do not understand this one completely
caller = dbstack;
if numel(caller)>=3 % If called from another script
    fprintf('**************** q-resolved ************************\n')
    fprintf('Adjusting data to reconstruct q-resolved index %05d\n',indsq)
    q_indices = projection(ii).par.qresolved_q{indsq};
    s.q_range_A = projection(1).integ.q(q_indices(:));
    p.add_name = [p.add_name sprintf('_qres_ind_%05d',indsq)];
    for ii = 1:numel(projection)
        projection(ii).data = projection(ii).data_q(:,:,:,indsq);
        % Cleaning up a bit, making sure that parameters that would be
        % incorrect cannot be used.
%         projection(ii).par.which_qrange = [];
%         projection(ii).par.r_sum = {};
    end
    doing_q_resolved = true;
    % Check for Nan
    for ii = 1:numel(projection)
        if any(~isfinite(projection(ii).data(:)))
            fprintf('Found NaN or non-finite value in projection %d\n',ii)
            figure(4); imagesc(any(~isfinite(projection(ii).data),3)); title('Found NaN')
            projection(ii).window_mask(any(~isfinite(projection(ii).data),3)) = 0; %#ok<SAGROW>
            projection(ii).data(~isfinite(projection(ii).data)) = 0; %#ok<SAGROW>
            figure(3); imagesc(projection(ii).window_mask); title(sprintf('Projection %d, Corrected mask to ignore pixel with non-finite value',ii))
%             return                                               
%             keyboard
        end
    end
%     keyboard
else
    % Adding q-range used in the integration to the tensor
    fprintf('Adjusting data \n')
    whichqrange = projection(1).par.which_qrange;
    q_indices = projection(1).par.r_sum{1,whichqrange};
    %s.q_range_A = projection(1).integ.q(q_indices(:)); % I think it will be useful to not have only the extreme values but also the number of q-indices used, so I keep all
    % IGNORED THIS VARIABLE AS IT IS NOT BEING USED. WEIRD...
    doing_q_resolved = false;
end

%% Run settings AD

p.mode = 1;                      % RSD: 1 for AD, 0 for symbolic
p.method = "bilinear";          % RSD: Choose method of interpolation.
p.filter_2D = 1;                % RSD: The best filter.
p.GPU = 0;
p.python = 1; %RSD: Improve safety here.
p.batch = 0; %RSD: ALWAYS

% Defining filenames for results. RSD: Moved down one cell to create
% add_name based on mode
if p.mode && p.GPU
    p.add_name = "AD_GPU";
elseif p.mode && p.python
    p.add_name = "AD_python";
elseif p.mode
    p.add_name = "AD";
else
    p.add_name = "symbolic";
end


make_3Dmask = 0;


%% Step 2.2 optimization of coefficients over the symmetric intensity: only a0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.opt_coeff = [1];             % parameters for the optimization

p.find_orientation = 0;        % Optimize over the SH angles: theta and phi (true or false)

p.regularization = 0;          % Sieves regularization on the coefficients. (true or false)
p.regularization_angle = 0;    % Regularization of the angles.  (true or false)
% =======
%optimize

%parameters for the sieves regularization (blurring the gradient)
%for the coefficients only
p.slice = 0; %in case the optimization is done on a slice (so it is using a 2D kernel)
kernel3D = window3(5,5,5,@hamming);
p.kernel = kernel3D./sum(kernel3D(:)); % for normalization (sum equals 1)

p.itmax = 10; %20                % maximum number of iterations: about 20
p.skip_projections = 1;         % = 1, for not skipping projections


sym_opt_filename = fullfile(p.optimization_output_path, sprintf('result_%s%s_symSAXS.mat', p.sample_name, p.add_name)) ;  %sprintf(fullfile('%s','result_%s%s_symSAXS.mat'), p.optimization_output_path, p.sample_name, p.add_name); 
%FIXED PATH: A GOOD RULE OF THUMB SHOULD BE TO HAVE FULLFILE AS THE OUTMOST FUNCTION WHEN DEFINING A PATH.
fprintf('filename: %s', sym_opt_filename);

p.avoid_wrapping=0;            % avoid wrapping over 2Pi

p.save.output_filename = sym_opt_filename;
p.save.image_filename = sprintf('%s/optimization_sym_int_%s', p.figures, p.add_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%optimize
fprintf('****** Step 2.2 optimization of coefficients over the symmetric intensity: only a0 ******\n')
tic
[p, s] = optimization.optimize_SH(projection, p, s); % RSD : true for AD, false for symbolic
toc
fprintf('Saving results in %s\n',p.save.output_filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2.3: Define a 3D mask (voxels that aren't going to be optimized)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % 1 to create a 3D mask, 0 to see the symmetric intensity tomogram
mask.cut_off = []; % = [] for automatic detection of the mean value using Otsu threshold
mask.diam_cylinder = []; % Diameter of a cylinder outside, = [] does not create the mask
mask.gauss_filter = 1; % Apply a filter to the data to reduce noise
mask.openning = 1;     % Apply an openning operation on the mask, this makes it 
                       % larger as is important to not have too tight support constraint.
coloraxis = 'auto'; % apply same color scale to the three projections xy, yz and xz fron view3
viewpar.interactive = 0;  % = true if you want to explore the volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Step 2.3: Define a 3D mask (voxels that aren''t going to be optimized)\n')

% load the result from symmetric intensity
fprintf('Loading %s\n',sym_opt_filename)
load(sym_opt_filename);

% needed for plotting
    
if make_3Dmask
   %[s.mask3D, cut_off] = create_mask3D(s.a(1).data, mask);
   s.mask3D = fasit.mask3D;
else
    s.mask3D = ones(size(s.a(1).data));
    cut_off = 0;
end

%save the updated 3Dmask in symmetric intensity
fprintf('Saving the updated 3D mask in %s\n', sym_opt_filename)
save(sym_opt_filename, 's', 'p', 'e', '-v6');

% plot the masked symmetriuc intensity reconstruction

%if numel(caller)<3 % If not called from another script
%    view3(s.mask3D.*s.a(1).data + max(s.a(1).data(:))*(1-s.mask3D), coloraxis, viewpar)
%    colormap(plotting.franzmap);
%end



%% Step 2.4: optimization of SH angles (phi and theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Step 2.4: optimization of SH angles (phi and theta)\n')
sym_int = importdata(sym_opt_filename);
fprintf('Loading %s\n',sym_opt_filename)

p = sym_int.p;

p.opt_coeff = [0, 0, 0]; 

p.find_orientation = 1;           % Optimize over the main orientation of the structure (true or false)

p.regularization = 0;             % Sieves regularization on the coefficients. (true or false)
p.regularization_angle = 0;       % Regularization of the angles. (true or false)
p.regularization_angle_coeff = 0; % mu for regularization of angle, needs to be found with L-curve

p.itmax = 20; %50; %30;           % maximum number of iterations: about 50
p.skip_projections = 1; % = 1, for not skipping projections


p.avoid_wrapping = 1;   % Avoid wrapping of the angle (to keep angle between 0 and 2pi) (true or false)

angle_opt_filename = sprintf('%s/result_%s_q%d-%d_angles_%s.mat', p.optimization_output_path, ...
    p.sample_name, projection(1).par.r_sum{1}(1),  projection(1).par.r_sum{1}(end), p.add_name);

p.save.output_filename = angle_opt_filename;
p.save.image_filename = sprintf('%s/optimization_angles_%s', p.figures, p.add_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do optimization
close all

% parameters for angles only
% values for initial guess of coeffiecient, relative to symmetric intensity
% coefficients of the SH
% Vector with spherical harmonic parameters
l = [0 2 4];  % Polar order
m = [0 0 0];  % Azimuthal order

%%%define ratio of coefficient (fixed in this step) for bone: 1, -3, 6,
%%%%RSD: Ratios should be changed. Have carbon knot.
a_ratio = [1, 5, 10];

for ii = 1:numel(l)
    s.a(ii).data = sym_int.s.a(1).data / a_ratio(ii);
    s.a(ii).l = l(ii);
    s.a(ii).m = m(ii);
end

%optimize
tic
[p, s] = optimization.optimize_SH(projection, p, s); %RSD: remember AD (p.mode)
toc


%RSD: Also, why does the error increase from a0? Perhaps because a2 and a4 are introduced by a guess. 

%% Step 2.5: optimization of SH coefficients: non symmetric
% optimize the other coefficients a2, a4 and a6, keeping a0 constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.opt_coeff = [0, 1, 1, 1];

p.find_orientation = 0;     % Optimize over the main orientation of the structure (true or false)
p.regularization = 0; %1;       % apply Sieves regularization on the coefficient (true or false)
p.regularization_angle = 0; % Regularization of the angles. (true or false)

%parameters for the sieves regularization (blurring the gradient)
p.slice = 0;
kernel3D = window3(5,5,5,@hamming);    % strength of Sieves regularization
p.kernel = kernel3D./sum(kernel3D(:)); % for normalization (sum equals 1)

p.itmax = 20; %20;            % maximum number of iterations: about 20
p.skip_projections = 1;  % = 1, for not skipping projections

p.avoid_wrapping = 0;    % avoid wrapping over 2Pi of the angle RSD: What should be chosen?

coef_a1const_filename = sprintf('%s/result_%s_q%d-%d_coef_a1const_%s.mat',p.optimization_output_path, ...
    p.sample_name, projection(1).par.r_sum{1}(1),  projection(1).par.r_sum{1}(end), p.add_name);

p.save.output_filename = coef_a1const_filename;
p.save.image_filename = sprintf('%s/optimization_coef_%s', p.figures, p.add_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do optimization
close all

fprintf('Step 2.5: optimization of SH coefficients: non symmetric\n')

fprintf('Loading symmetric intensity results from %s\n',sym_opt_filename)
sym_opt = importdata(sym_opt_filename);

% parameters for optimization of coefficients with constant a1 fixed from initial optimization
% loas parameters from the angle
fprintf('Loading parameters from %s\n',angle_opt_filename)
angle_opt = importdata(angle_opt_filename); % RSD: Cannot see if the optimised angles are loaded. But believe they remain in memory when everything is run at once. TODO: Make another run file. 

p.phi_det = angle_opt.p.phi_det;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vector with spherical harmonic parameters. Setting an initial guess for
% the coefficients.
l = [0 2 4 6];  % Polar order
m = [0 0 0 0];  % Azimuthal order
a = [1 0.2 0.1 0.05];   % Coefficients

for ii = 1:numel(l)
    if p.opt_coeff(ii)
        s.a(ii).data = ones(p.ny,p.nx,p.nz)*a(ii);
    else      % This is only needed in case of loading from the file. I am not sure I like this idea that we always start each step by loading from file by default.
        s.a(ii).data = sym_opt.s.a(1).data;
    end
    s.a(ii).l = l(ii);
    s.a(ii).m = m(ii);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%optimize
tic
[p, s] = optimization.optimize_SH(projection, p, s);
toc
%% Step 2.6: find optimal regularization parameters for final step
%RSD: Not investigated at all yet. 

%% Step 2.7: final optimization: combine all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine an approximation of the regularization coefficient
load(coef_a1const_filename)

t = sum(((s.a(1).data.*optimization.spherical_harm(0,0,0,0)).^2),4); %plot the first coefficent
tt = (sum(t, 3).*projection(1).window_mask); %axis tight equal xy
mu_guess = mean(mean(tt));

% optimize all coefficients (a0, a2, a4 and a6) and angles (theta and phi)
% whihout any constriction
p.opt_coeff = [1, 1, 1, 1];

p.find_orientation = 1;      % Optimize over the main orientation of the structure

%RSD: Ignore regularization for now.
p.regularization = 0;%1;        % Sieves regularization on coeff
p.regularization_angle = 0;% 1; %regularization of angle
p.regularization_angle_coeff = mu_guess;  %find the appropriate mu with L-curve from step 2.6 % RSD: Is this correct?

%parameters for the sieves regularization (blurring the gradient)
kernel3D=window3(5,5,5,@hamming);
p.kernel=kernel3D./sum(kernel3D(:)); % for normalization (sum equals 1)

p.itmax = 50; %50;                          % maximum number of iterations: about 50-100
p.skip_projections = 1;                % = 1, for not skipping projections

p.avoid_wrapping = 1;     % avoid wrapping over 2Pi of the angle

all_again_filename = sprintf('%s/result_%s_q%d-%d_all_again_%s.mat', p.results, ...
    p.sample_name, projection(1).par.r_sum{1}(1),  projection(1).par.r_sum{1}(end), p.add_name);

p.save.output_filename = all_again_filename;
p.save.image_filename = sprintf('%s/optimization_all_q%d-%d_%s', p.figures, projection(1).par.r_sum{1}(1),  projection(1).par.r_sum{1}(end), p.add_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

fprintf('Step 2.7: final optimization: combine all\n')

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
tic
[p, s] = optimization.optimize_SH(projection, p, s);
toc

