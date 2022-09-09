% This template is used to create fourier analysis images, and
% SASTT loading data analysis
close all
clear variables
addpath ../base/
%% STEP 1: give the locations to find the raw data and to save the processed data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.online_analysis = 1;            %  = 1 for online analysis, =0 for offline, needs user_name
par.sample_name = 'sample_name';    %  name for the subfolders
user_name = '17041'; % e-account number, without 'e', used only if par.online_analysis = 0
fmt.parallel_plotting = false;      % = true will plot using a parfor, 
                                    % plots will not be shown on the screen
par.use_nexus = 1;               % If NeXus raw data is available then all the information
                                    % can be obtained from it, such as mcs,
                                    % sgalil positions, spec data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if par.online_analysis
    % use online (e-account)
    par.user_name = beamline.identify_eaccount;
    base_dir_beamtime = '~/Data10/';
    base_dir_work = base_dir_beamtime;
else % for off-line analysis: Daas
    par.user_name = ['e' user_name];
    % use offline (p-group)
    base_dir_beamtime = sprintf('/das/work/p%s/p%s/online_data', user_name(1:2), user_name);
    base_dir_work = sprintf('/das/work/p%s/p%s/', user_name(1:2), user_name);
end
fnames.data_dir  = sprintf('%s/analysis/', base_dir_beamtime); % Where the radial integration data is located. Read only
fnames.save_dir  = sprintf('%s/analysis/', base_dir_work); % where to save the processed data

if par.use_nexus % if we use Nexus the parameters below are used as filter in io.nexus_read
    fnames.mcs_dir   = 'mcs';       % which group for mcs data, read only
    fnames.pos_file  = 'sgalil';    % which group for positions data, 'sgalil' or 'spec' to get positions from spec data or spec header
    fnames.spec_file = 'specES1';   % which group for spec data
    fnames.base_dir_beamtime = base_dir_beamtime;
else
    fnames.mcs_dir   = sprintf('%s/mcs/', base_dir_beamtime); % where to look for the mcs data, read only
    fnames.pos_file  = sprintf('%s/sgalil/S%%05d.dat', base_dir_beamtime); % where to look for the positions data, if using owis you can get positions from spec by putting here the SPEC dat-file. Read only
    fnames.spec_file = sprintf('%s/', base_dir_beamtime);
end
% directory-name part used relative to basedir_fourier, should include an ending slash
fnames.sub_dir = [par.sample_name '/'];
% Load the data
% directory for loading integrated data, empty to choose the default
fnames.basedir_integrated_data = [ fnames.data_dir 'radial_integration/'  ]; % might need to readjust in case /S10000-19999/
fnames.fext = '.h5';
fnames.basename_integrated_data = [par.user_name '_1_%05d_00000_00000_integ' fnames.fext];
% save Fourier components to this file, empty to choose the default
fnames.basedir_fourier = [];
fnames.filename_fourier = [];
% save figures to this directory, empty to choose the default
fnames.basedir_figures = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STEP 2: perform background correction?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fnames.bgr_scan_no = [ ]; % =[] for no correction, otherwise background scan number, example: 1120
% the following options are relevant only, if bgr_scan_no is not empty:
fnames.bgr_point_range = [ ]; % = average over all points of the specified scan, or select range, example: 1:201;
fnames.basename_bgr_data = [ ]; % directory in which the background scan data are located and file-name template to insert the scan number into. E.g. [ fnames.data_dir 'radial_integration_SAXS/e' user_name '_1_%05d_00000_00000_integ.mat' ];
% outlyer background points (dust, scratch, etc.) can be filtered out using a simple median -filtering
fnames.bgr_max_dev_from_median = 0.10; % =[] to disable median filtering of outlyers, otherwise discard points with deviation from median value higher than this factor, example: 0.10 for max. 10% deviation
fnames.bgr_median_radii_range = 20:200; % =[] for sum over all radii, not relevant if bgr_max_dev_from_median is empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STEP 3: perform transmission correction?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transmisison correction data recorded using the MCS,
% (four_dat.f1_amp and four_dat.f2_amp will be corrected for the transmission data,
% and the transmission data will be stored on four_dat.trans_data)
% transmission data are taken from this channel (column) of the MCS data
fnames.trans_mcs_channel = 3;
fnames.basedir_trans_data = [fnames.data_dir 'online/stxm/data/']; %=[] for no correction
% first and last scan no. are set into this name
fnames.basename_trans_data = 'stxm_scans_%05d-%05d_mcs.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% This section is used to see a scattering curve to determine r_sum for next step
find_rsum_values = 0; % = 1 to show integrated pattern to select ranges, = 0 to ignore
if find_rsum_values == 1
    par.test_scan = 48;       % which scan number is used to check r_sum
    par.point_range = [];     % = [ ] to average over all points in the scan line/loopscan or = [a:b] select a range of points
    par.PlotMulQ = 1;         % = 0, don't multiply. = 1 multiply I(q) by a q to the power 
    par.QMulPow = 0.5;        % multiply intensity with q to the power of this value, default is [ ] for no multiplication
    [par] = find_rsum(par, fnames);
end

%% STEP 4: define the q-ranges of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define which radii ranges shall be integrated (q-resolved)
% choose the desired q-range to perform SASTT in step 8 by selecting par.which_qrange = ?;
par.r_sum = {19:24, 57:68, 82:378};

fnames.save_name_suffix = {'first_range', 'second_range', 'third_range' };
if (length(fnames.save_name_suffix) ~= length(par.r_sum))
    error('The number of elements in fnames.save_name_suffix is inconsistent with the number of summation ranges in par.r_sum');
end

% Background points from to, 0 0 for no background subtraction.
% with an r_sum area from 50:60 and r_bgr 5 2 the background is averaged
% from (50-5) to (50-2) and from (60+2) to (60+5)
par.r_bgr = [ 0 0; 0 0; 0 0 ];
if (size(par.r_bgr,1) ~= length(par.r_sum))
    error('The number of elements in par.r_bgr is inconsistent with the number of summation ranges in par.r_sum');
end

% power law exponent, only relevant for background interpolation for
% subtraction, -1 means determination by using the intensities in the above
% specified background region, otherwise 0 to 4 as physical meaningful,
% hard-coded values
par.I_bgr_pow = [ -1, -1, -1 ];
if (length(par.I_bgr_pow) ~= length(par.r_sum))
    error('The number of elements in par.I_bgr_pow is inconsistent with the number of summation ranges in par.r_sum');
end

% if non-zero and background area specified via par.r_bgr, 
% then try peak analysis:
% 1 - simple peak analysis (try this for online analysis)
% 2 - fitting a Gaussian peak (slower than the simple analysis)
par.peakfit = [ 0, 0, 0 ];
if (length(par.peakfit) ~= length(par.r_sum))
    error('The number of elements in par.peakfit is inconsistent with the number of summation ranges in par.r_sum');
end

% provide these parameters in case q is not saved in radial integration
% parameters of the setup
% par.pix_size = 0.172; %pixel size
% par.det_dist = 2.1604e3;  % detector distance
% par.lambda_A = 12.398/11.4; % 12.398/Energy (keV)

%% STEP 5: Tasks for the code and plotted output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define which steps have to be done
processing.load_data = 0; % 1 to force re-load and re-process
processing.save_fig  = 1; % 1 to save the images
processing.plot_data = 1; % 1 to avoid closing and opening the windows
processing.movie_fig = 0; % 1 to save the images as a video
processing.print_fig = 0; % 1 to print the figures (beamline printer as default)
% choose which figures to plot
fmt.fig_sym_int = 1;        % symmetric intensity over all segments: f1_amp
fmt.fig_asym_int = 1;       % amplitude of the anisotropy: f2_amp
fmt.fig_orientation = 1;    % angle of orientation (in degrees): f2_phase
fmt.fig_orient_degree = 1;  % degree of orientation from f2_amp/f1_amp
fmt.fig_orient_degree_histogram = 0; % histogram of the degree of orientation vs. orientation
fmt.fig_asym_sym_int = 1;   % combination of orientation (colorwheel), degree of orientation (hue)
fmt.fig_scatt_asym = 0;     % scattering asymmetry
fmt.fig_cos_dev = 0;        % how the fit of azimuthal integration deviates from a cosine
fmt.fig_I_pow = 0;          % slopes of I(q) in the selected q-region
% simple peak analysis
fmt.fig_simple_peak.pos = 1;       % peak position
fmt.fig_simple_peak.ampl = 1;      % peak amplitude
fmt.fig_simple_peak.width = 1;     % peak width
fmt.fig_simple_peak.all = 1;       % peak amplitude, position, width
% peak fitting
fmt.fig_fit_peak.pos = 1;       % peak position
fmt.fig_fit_peak.ampl = 1;      % peak amplitude
fmt.fig_fit_peak.width = 1;     % peak width
fmt.fig_fit_peak.all = 1;       % peak amplitude, position, width
% colormap scale: scale the caxis of each figure +/- this values times the standard
% deviation of the linear data
fmt.std_scale = 2.0;
% if the following fields are undefined or empty (i.e., [] ) then auto-scaling is
% applied to a maximum value calculated as the median value plus
% fmt.std_scale times the standard deviation of the values.
fmt.asymmetric_amplitude_scale_max = [];
fmt.symmetric_amplitude_scale_max = [];
% interpolation (smaller pixels) of the plotted data, 0 means off, greater than 1 will
% become very slow
fmt.interpol = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STEP 6: segmentation
par.prepare_segmentation = 0; % =0 to switch segmentation off, = 1 to perform the segmentation


%% STEP 7: SASTT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare data for SASTT
par.prepare_SASTT = 0;  % = 1 to prepare data to SASTT
par.which_qrange = 3;   % which q_range from par.r_sum to load for SASTT optimization
par.tomo_axis_x = 0;    % = 1 if the tomographic axis is horizontal
    % prepare data for q_resolved SASTT analysis
    perform_qres = 0; %=1 to perform q resolved, = 0 to ignore next subsection
    par.qresolved_q = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if perform_qres == 1 && par.prepare_SASTT == 1
        par.pixel_range = 5:650; % range of pixels usually from 5:900 (exclude noisy values)
        par.test_scan = 28;      % which scan number to visualized
        par.point_range = [];    % which point range to consider (range of points in the line): = [] for all
        par.int_integ = 40;      % how many q intervals should be considered for SASTT. The number is approximate, 
                                 %    if logscale is selected along with many intervals then bins that are repeated in low q are removed.
        par.log_scale = 0;       % = 0 for sampling q linearly and = 1 for sampling q in log scale
        plot_qres = 0;           % = 0: if not plotting is requested, = 1: if plotting is requested   
        [par] = find_qrange(par, fnames, plot_qres);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STEP 8: Choose the data to analyze

% frames = get_frameID(base_dir_beamtime, samplename','S46290');
% frames = get_frameID(base_dir_beamtime, 'Dataset_id',45);

% Or you can simply analyze based on known frameid numbers:
frames = [42:252];
wait_for_next_frame = false;

%
for ii=1:numel(frames)
    if wait_for_next_frame
        waitnow = true;
        while waitnow
            [scann_next,p_next] = get_scannum_frame(base_dir_beamtime,frames(ii)+1);
            if ~isempty(scann_next)
                waitnow = false;
            else
                fprintf('Waiting for next frame, pausing 2 seconds\n')
                pause(2)
            end
        end
    end
    [scann,p] = get_scannum_frame(base_dir_beamtime,frames(ii));
    fnames.rel_dir = '';
    % number of the first scan for this sample
    fnames.first_scan_no = scann(1);
    fnames.frame_id = p.frame_id;
    % scan dimension in intervals (not points)
    par.ind_max_x = p.saxs_intervals(1);
    par.ind_max_y = p.saxs_intervals(2);
    fnames.scan_no_step = 1;
    % pixel size in mm
    par.x_scale = p.saxs_stepsize(1);
    par.y_scale = p.saxs_stepsize(2);
    % the horizontal x-axis is the fast scan axis
    par.fast_axis_x = 1-p.fast_axis_y;
    par.snake_scan = p.snake_scan;
    % perform the analysis hard coded in this sub-routine
    fnames.rel_dir = sprintf('%s_scann%d',p.samplename, scann(1));
    analyze_one_scan(fnames,par,fmt,processing);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Legacy for scans before Oct 2019 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For used the old sSAXS scanning macro you dont have defined a "frame"
% in the metadata, in that case the following section should be copied from 
% spec_scanning_SAXS_caller_template_yyyy_mm_dd__hh_mm_ss.m
% This file is written by spec when the scanning_saxs macro is used
% It's saved in ~/Data10/cxs_software

% The format looks something like this

% SAXS data
% %% relative directory identifying the sample
% fnames.rel_dir = 'sample_name/';
% % number of the first scan for this sample
% fnames.first_scan_no = 471;
% % scan dimension in intervals (not points)
% par.ind_max_x = 100;
% par.ind_max_y = 200;
% fnames.scan_no_step = 1;
% % pixel size in mm
% par.x_scale = 0.040;
% par.y_scale = 0.040;
% % the horizontal x-axis is the fast scan axis
% par.fast_axis_x = 0;
% par.snake_scan = 1;
% % perform the analysis hard coded in this sub-routine
% analyze_one_scan(fnames,par,fmt,processing);




%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
%
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing la this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS scanning SAXS package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% Additionally, any publication using the package, or any translation of the 
%     code into another computing language should cite:
%    O. Bunk, M. Bech, T. H. Jensen, R. Feidenhans'l, T. Binderup, A. Menzel 
%    and F Pfeiffer, “Multimodal x-ray scatter imaging,” New J. Phys. 11,
%    123016 (2009). (doi: 10.1088/1367-2630/11/12/123016)
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
