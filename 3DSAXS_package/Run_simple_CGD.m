
%% Step 2.1a Load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent = cd;
base_path = '/Data sets/';  % = '~/Data10/' for online analysis,
base_path = [parent base_path] ;
                                            % provide the path for offline analysis
                                            % Ex: '/das/work/p16/p16649/'
sample_name = 'Synthetic_sample_ny4_ny4_all_coeffs'; % 'SASTT_carbon_knot_aligned_ASTRA_correctedny4nx4'; %'Synthetic_sample_ny4_ny4_all_coeffs'; %'sample_name';     % name given in the saxs_caller_template
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
p.numsegments = length(projection(1).integ.phi_det)/2;

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

l = [0 2 4 6];  % Polar order
m = [0 0 0 0];  % Azimuthal order
base = 0.69;
a = [base  base*1/5 base*1/48 base*1/96];   % Coefficients

%%% Initial values for the SH angles
theta_init = pi/2; %45 degrees:  Polar orientation of structure
phi_init = pi/2; %45 degrees: Azimuthal orientation of structure

% calculate the initial values
s.theta.data = ones(p.ny,p.nx,p.nz)*theta_init;  % Polar orientation of structure
s.phi.data   = ones(p.ny,p.nx,p.nz)*phi_init;      % Azimuthal orientation of structure

% Choose orientation pattern
theta_op = ones(p.ny,p.nx,p.nz) ;
phi_op = ones(p.ny,p.nx,p.nz) ;

% This model is domains
theta_op(1: p.ny/2 , 1: p.nx/2, 1: p.nz/2 ) = 0;
theta_op(1:p.ny/2 , 1 + p.nx/2: p.nx, 1: p.nz/2) = pi/4;
theta_op(1+p.ny/2: p.ny , 1 : p.nx/2, 1: p.nz/2) = pi/3;
theta_op(1+p.ny/2: p.ny , 1 + p.nx/2: p.nx, 1: p.nz/2) = pi/2;

theta_op(1: p.ny/2 , 1: p.nx/2, 1+ p.nz/2 : p.nz) = 0 + pi/2;
theta_op(1: p.ny/2 , 1 + p.nx/2: p.nx, 1+ p.nz/2 : p.nz) = pi/4 + pi/2;
theta_op(1+ p.ny/2: p.ny , 1 : p.nx/2, 1+ p.nz/2 : p.nz) = pi/3 + pi/2;
theta_op(1+ p.ny/2: p.ny , 1 + p.nx/2: p.nx,1+ p.nz/2 : p.nz) = pi;

phi_op(1: p.ny/2 , 1: p.nx/2, 1: p.nz/2) = 0;
phi_op(1: p.ny/2 , 1 + p.nx/2: p.nx, 1: p.nz/2) = pi/4;
phi_op(1+p.ny/2: p.ny , 1 : p.nx/2, 1: p.nz/2) = pi/3;
phi_op(1+p.ny/2: p.ny , 1 + p.nx/2: p.nx, 1: p.nz/2) = pi/2;

phi_op(1: p.ny/2 , 1: p.nx/2, 1+ p.nz/2 : p.nz) = 0 + pi/2;
phi_op(1: p.ny/2 , 1 + p.nx/2: p.nx, 1+ p.nz/2 : p.nz) = pi/4 + pi/2;
phi_op(1+p.ny/2: p.ny , 1 : p.nx/2, 1+ p.nz/2 : p.nz) = pi/3 + pi/2;
phi_op(1+p.ny/2: p.ny , 1 + p.nx/2: p.nx,1+ p.nz/2 : p.nz) = pi;

s.theta.data = theta_op;  % Polar orientation of structure
s.phi.data   = phi_op;      % Azimuthal orientation of structure

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
p.filter_2D = 3;                % RSD: The best filter.
p.GPU = 0;
p.python = 0; %RSD: Improve safety here.

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


% needed for plotting
    
if make_3Dmask
   [s.mask3D, cut_off] = create_mask3D(s.a(1).data, mask);
else
    s.mask3D = ones(size(s.a(1).data));
    cut_off = 0;
end

%save the updated 3Dmask in symmetric intensity
% fprintf('Saving the updated 3D mask in %s\n', sym_opt_filename)
% save(sym_opt_filename, 's', 'p', 'e', '-v6');



%% Step 2.7: final optimization: combine all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine an approximation of the regularization coefficient
%load(coef_a1const_filename)

% optimize all coefficients (a0, a2, a4 and a6) and angles (theta and phi)
% whihout any constriction
p.opt_coeff = [1, 1, 1, 1];

p.find_orientation = 1;      % Optimize over the main orientation of the structure

%RSD: Ignore regularization for now.
p.regularization = 0;%1;        % Sieves regularization on coeff
p.regularization_angle = 0;% 1; %regularization of angle

%parameters for the sieves regularization (blurring the gradient)
kernel3D=window3(5,5,5,@hamming);
p.kernel=kernel3D./sum(kernel3D(:)); % for normalization (sum equals 1)

p.itmax = 10; %50;                          % maximum number of iterations: about 50-100
p.skip_projections = 1;                % = 1, for not skipping projections

p.avoid_wrapping = 1;     % avoid wrapping over 2Pi of the angle

all_again_filename = sprintf('%s/result_%s_q%d-%d_CGD_%s.mat', p.results, ...
    p.sample_name, projection(1).par.r_sum{1}(1),  projection(1).par.r_sum{1}(end), p.add_name);

p.save.output_filename = all_again_filename;
p.save.image_filename = sprintf('%s/optimization_all_q%d-%d_%s', p.figures, projection(1).par.r_sum{1}(1),  projection(1).par.r_sum{1}(end), p.add_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

fprintf('Step 2.7: final optimization: combine all\n')

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

