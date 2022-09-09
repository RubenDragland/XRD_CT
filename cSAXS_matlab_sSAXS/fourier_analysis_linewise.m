% [ four_dat, data_process] = fourier_analysis_linewise(fnames,par,fmt)

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

function [ four_dat, data_process] = fourier_analysis_linewise(fnames,par,fmt)
import io.mcs_mesh
import utils.adjust_projection
import utils.fopen_until_exists

% in case its not specified, it wil consider that snake_scan = 0
if ~isfield(par,'snake_scan')
    par.snake_scan = 0;
end
if ~isfield(par,'prepare_segmentation')
    par.prepare_segmentation = 0;
end
if ~isfield(par,'prepare_SASTT')
    par.prepare_SASTT = 0;
end
if ~isfield(par,'which_qrange')
    par.which_qrange = [];
end
if ~isfield(par,'fast_axis_x')
    par.fast_axis_x = 0;
end

fprintf('loading the integrated data\n');

if ((~isfield(fnames,'scan_no_step')) || (isempty(fnames.scan_no_step)))
    scan_no_step = 1;
else
    scan_no_step = fnames.scan_no_step;
end

if ~isfield(fnames, 'fext')
    fnames.fext = '.h5';
end

% the fast and slow axes depend on the orientation of the scan.
if (par.fast_axis_x)
    ind_max_fast = par.ind_max_x;
    ind_max_slow = par.ind_max_y;
else
    ind_max_fast = par.ind_max_y;
    ind_max_slow = par.ind_max_x;
end

% load the first file to get the ranges ind_max_x, ind_max_y
format_05 = -1;
[ filename, format_05 ] = ...
    filename_integrated_intensities(fnames.basedir_integrated_data, ...
    fnames.basename_integrated_data, ...
    fnames.first_scan_no,format_05, fnames.fext);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wait for the file to be available: important for online analysis
wait_ind = 0;
while (~exist(filename,'file'))
    if (wait_ind == 0)
        fprintf('Waiting for %s to appear.\n',filename);
    end
    pause(60);
    wait_ind = wait_ind +1;
end
if (wait_ind > 0)
    pause(60);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the data and check the range of the radii used
if strcmpi(fnames.fext,'.h5')
    first_data = io.HDF.hdf5_load(filename);
else 
    first_data = load(filename);
end
ind_r_max = length(par.r_sum);
num_segments = size(first_data.I_all,2);

% store radial-integration normation sums (no. of pixels),
% q-values, 
% and positions of angular segments
integ.norm_sum = first_data.norm_sum;
if isfield(first_data, 'q')
    integ.q = first_data.q;
else
    fprintf('Radial_integration DOES NOT contain q-values\n');
end
if isfield(first_data, 'phi_det')
    integ.phi_det = first_data.phi_det;
else
    fprintf('Radial_integration DOES NOT contain phi_det\n');
    % try backward compatibility -- orientation setting depends on this
    integ.phi_det = ((360/num_segments)/2:(360/num_segments):360);  % 11.25:22.5:348.75;
end

%%
%define the type of scan: mesh (one scan number) or linewise (one scan
%number per line along the fast axis)
if (size(first_data.I_all,3) >= (ind_max_fast+1)*(ind_max_slow+1))
    fprintf('Switching to mesh scan mode (one data file for the whole scan).\n');
    mesh_scan_mode = 1;
else
    mesh_scan_mode = 0;
end
%%
% initialize the output
f1_amp   = zeros(ind_max_fast+1,ind_max_slow+1,ind_r_max);
f2_amp   = f1_amp;
f2_phase = f1_amp;
I_dev    = f1_amp;
I_cos_dev= f1_amp;
I_pow    = f1_amp;
% fit peak profiles
fit_peak_pos   = f1_amp;
fit_peak_ampl  = f1_amp;
fit_peak_width = f1_amp;
% simple and faster peak analysis
simple_peak_pos   = f1_amp;
simple_peak_ampl  = f1_amp;
simple_peak_width = f1_amp;
% for SAS TT
data_process = {};
if (par.prepare_segmentation)
    data_process.sym_data = zeros(ind_max_fast+1,ind_max_slow+1,size(integ.norm_sum,1));
    if (~isempty(fnames.basedir_trans_data))
        data_process.trans_data = zeros(ind_max_fast+1,ind_max_slow+1);
    end
end
if (par.prepare_SASTT)
    data_process.sector_data = zeros(ind_max_fast+1,ind_max_slow+1,num_segments/2);
    integ.norm_sum_avg = zeros(num_segments/2,length(par.r_sum{par.which_qrange}));
    % temporary variable for parfor loop
    integ_norm_sum_avg = zeros(ind_max_slow+1,num_segments/2,length(par.r_sum{par.which_qrange}));
    if (isfield(par,'qresolved_q') && (~isempty(par.qresolved_q)))
        data_process.sector_data_q = zeros(ind_max_fast+1,ind_max_slow+1,num_segments/2,numel(par.qresolved_q));
        integ.norm_sum_avg_q = zeros(num_segments/2,max(cellfun(@numel,par.qresolved_q)),numel(par.qresolved_q));
        % temporary variable for parfor loop
        integ_norm_sum_avg_q = zeros(ind_max_slow+1,num_segments/2,max(cellfun(@numel,par.qresolved_q)),numel(par.qresolved_q));
    end
end

% do background subtraction: average and not point by point
if (~isempty(fnames.bgr_scan_no))
    [ filename_bgr, format_05 ] = ...
        filename_integrated_intensities('',fnames.basename_bgr_data, ...
        fnames.bgr_scan_no,format_05, fnames.fext);
    % load the background
    fprintf('Loading background data from %s.\n',filename_bgr);
    bgr_data = load(filename_bgr, 'I_all', 'norm_sum'); % load only varaibles that are used
    % determine over the specified point range
    % averaged background intensity,
    % keep unavailable intensities flagged as -1
    title_str = strrep(filename_bgr,'\','\\');
    title_str = strrep(title_str,'_','\_');
    I_bgr_all = calc_I_point_avg(bgr_data.I_all,fnames,title_str);
else
    % no background correction
    I_bgr_all = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correct with transmission data:
%
data_adjusted =[];
if (~isempty(fnames.basedir_trans_data))||(~isempty(fnames.pos_file))
    if isempty(fnames.basedir_trans_data)
        positions_only = 1;
    else
        positions_only = 0;
    end
    % load mcs_data
    % allow for MCS data being recorded in a separate scan of identical
    % dimensions by optional use of an offset for the file number (example:
    % WAXS data without beam-stop diode are corrected for transmission
    % recorded in SAXS, assuming 1:1 correspondence of the sample
    % positions)
    if (isfield(fnames,'trans_mcs_file_no_offs'))
        trans_mcs_file_no_offs = fnames.trans_mcs_file_no_offs;
    else
        trans_mcs_file_no_offs = 0;
    end
    
    if par.use_nexus
        [mcs_data, data_adjusted, pos_data] = mcs_mesh(fnames.first_scan_no+trans_mcs_file_no_offs, ...
            ind_max_slow,'ChToPlot',fnames.trans_mcs_channel, ...
            'FastAxisX', par.fast_axis_x, 'SnakeScan', par.snake_scan, 'UseNexus', par.use_nexus,...
            'FigureDir', [fnames.save_dir 'online/stxm/figures/'], 'DataDir', [fnames.save_dir 'online/stxm/data/'],...
            'DirBase', fnames.base_dir_beamtime ,'NexusFilterMCS', fnames.mcs_dir,  'NexusFilterPositions', fnames.pos_file, 'Positions_only', positions_only);
    else
        [mcs_data, data_adjusted, pos_data] = mcs_mesh(fnames.first_scan_no+trans_mcs_file_no_offs, ...
            ind_max_slow,'ChToPlot',fnames.trans_mcs_channel, ...
            'FastAxisX', par.fast_axis_x, 'SnakeScan', par.snake_scan, 'UseNexus', par.use_nexus, 'FnameBase', [par.user_name '_'],...
            'FigureDir', [fnames.save_dir 'online/stxm/figures/'], 'DataDir', [fnames.save_dir 'online/stxm/data/'],...
            'DirBase', fnames.mcs_dir,  'Pos_file', fnames.pos_file, 'Positions_only', positions_only);
    end
    trans_data = data_adjusted.transm;
    
    if size(f1_amp,1) < size(trans_data,1)
        warning('MCS array size is larger than expected, removing one point per line')
        trans_data = trans_data(1:end-1,:,:);
        
        data_adjusted.transm        = data_adjusted.transm(1:end-1,:,:);
        data_adjusted.positions_out = data_adjusted.positions_out(1:end-1,:,:);
        data_adjusted.scan_num      = data_adjusted.scan_num(1:end-1,:,:);
        if par.snake_scan
            aux         =  zeros(size(trans_data));
            aux_pos     =  zeros([size(aux) 2]);
            if find(data_adjusted.scan_point(:,1,1) == min(data_adjusted.scan_point(:,1,1))) == 1
                aux(:,1:2:end,:)        = data_adjusted.scan_point(1:end-1,1:2:end,:);
                aux(:,2:2:end,:)        = data_adjusted.scan_point(2:end,2:2:end,:);
                aux_pos(:,1:2:end,:)    = pos_data(1:end-1,1:2:end,:);
                aux_pos(:,2:2:end,:)    = pos_data(2:end,2:2:end,:);
            else
                aux(:,1:2:end,:) = data_adjusted.scan_point(2:end,1:2:end,:);
                aux(:,2:2:end,:) = data_adjusted.scan_point(1:end-1,2:2:end,:);
                aux_pos(:,1:2:end,:)    = pos_data(2:end,1:2:end,:);
                aux_pos(:,2:2:end,:)    = pos_data(1:end-1,2:2:end,:);
            end
            data_adjusted.scan_point = aux;
            pos_data = aux_pos;
        else
            pos_data = pos_data(1:end-1,:,:);
            data_adjusted.scan_point    = data_adjusted.scan_point(1:end-1,:,:);
        end
    end
    
    % histogram approach to determine the maximum intensity of the mcs file:
    % this is used for the relative transmission
    % use one million intensity bins over the intensity range
    int_min = min(min(trans_data));
    int_max = max(max(trans_data));
    int_step = (int_max-int_min) / 2^20;
    % calculate the bin indices from the intensities
    hist_ind = round((trans_data(:)-int_min) / int_step) +1;
    % increase each bin an index is pointing to by one
    hist_bins = zeros(1,max(hist_ind)-min(hist_ind)+1);
    for (ind_hist = 1:numel(hist_ind))
        hist_bins(hist_ind(ind_hist)) = hist_bins(hist_ind(ind_hist)) +1;
    end
    % calculate the cumulative sum
    hist_bins = cumsum(hist_bins);
    % find the indices to the bins with more pixels than the threshold
    n_high = find(hist_bins > 0.99 * numel(trans_data),1,'first');
    if (n_high > length(hist_bins))
        n_high = length(hist_bins);
    end
    norm_val = ((n_high-1) * int_step + int_min);
    fprintf('Normalizing channel %d by %.0f\n',fnames.trans_mcs_channel, ...
        norm_val);
    trans_data = trans_data ./ norm_val;
    % save in the results
    four_dat.trans_data = trans_data;
    % The original mcs data is used for the correction of SAXS patterns
    % individually
    % Here we should still check if a flip is needed
    if ~isempty(mcs_data)
        original_trans_data = squeeze(mcs_data(fnames.trans_mcs_channel,1,:,:))./ norm_val;
    else
        original_trans_data = [];
    end
else
    pos_data = [];
    original_trans_data = [];
end

% save the calculated normalized transmission value for every pixel
if (par.prepare_SASTT || par.prepare_segmentation)
    if (~isempty(fnames.basedir_trans_data))
        %save the air transmission calculated
        data_process.trans_air = single(norm_val);
    end
%     counter = 1;
end

% initialize temporary arrays
if (isfield(data_process,'sym_data'))
    data_process_sym_data = data_process.sym_data;
end
if (isfield(data_process,'trans_data'))
    data_process_trans_data = data_process.trans_data;
end
if (isfield(data_process,'sector_data_q'))
    data_process_sector_data_q = data_process.sector_data_q;
end
if (isfield(data_process,'sector_data'))
    data_process_sector_data = data_process.sector_data;
end


% parfor-loop, use for-loop instead to debug, plot intermediate results etc.
parfor (ind_slow=0:ind_max_slow) % go along each scan line
    if (rem(ind_slow,10) == 0)
        fprintf('%5d / %d\n',ind_slow+1,ind_max_slow+1);
    end
    
    % process one line of the raster scan
    [f1_amp_line, f2_amp_line, f2_phase_line, I_cos_dev_line, I_dev_line, I_pow_line, ...
        fit_peak_line, ...
        simple_peak_line, ...
        data_process_line] = ...
        process_one_line(par, fnames, mesh_scan_mode, ...
        fnames.first_scan_no + ind_slow*scan_no_step, format_05, ...
        ind_slow, ind_max_slow, ind_max_fast, ...
        I_bgr_all, integ, ...
        original_trans_data);
    
    % save return values in the corresponding arrays

    % Fourier analysis
    f1_amp(:,ind_slow+1,:) = f1_amp_line;
    f2_amp(:,ind_slow+1,:) = f2_amp_line;
    f2_phase(:,ind_slow+1,:) = f2_phase_line;
    % debug value
    I_cos_dev(:,ind_slow+1,:) = I_cos_dev_line;
    I_dev(:,ind_slow+1,:) = I_dev_line;
    % exp;onent of intensity decay
    I_pow(:,ind_slow+1,:) = I_pow_line;
    % peak analysis / fitting
    fit_peak_pos(:,ind_slow+1,:) = fit_peak_line.pos;
    fit_peak_ampl(:,ind_slow+1,:) = fit_peak_line.ampl;
    fit_peak_width(:,ind_slow+1,:) = fit_peak_line.width;
    simple_peak_pos(:,ind_slow+1,:) = simple_peak_line.pos;
    simple_peak_ampl(:,ind_slow+1,:) = simple_peak_line.ampl;
    simple_peak_width(:,ind_slow+1,:) = simple_peak_line.width;
    % SAS TT
    if (~isempty(data_process_line))
        if (isfield(data_process_line,'sym_data'))
            data_process_sym_data(:,ind_slow+1,:) = data_process_line.sym_data;
        end
        if (isfield(data_process_line,'trans_data'))
            data_process_trans_data(:,ind_slow+1) = data_process_line.trans_data;
        end
        if (isfield(data_process_line,'sector_data_q'))
            data_process_sector_data_q(:,ind_slow+1,:,:) = data_process_line.sector_data_q;
        end
        if (isfield(data_process_line,'sector_data'))
            data_process_sector_data(:,ind_slow+1,:) = data_process_line.sector_data;
        end
        if (isfield(data_process_line,'integ_norm_sum_avg'))
            % These normation sums should not depend on the line-index
            % ind_slow, i.e. setting these variables should be taken out of
            % the loop.
            integ_norm_sum_avg(ind_slow+1,:,:) = data_process_line.integ_norm_sum_avg;
        end
        if (isfield(data_process_line,'integ_norm_sum_avg_q'))
            % These normation sums should not depend on the line-index
            % ind_slow, i.e. setting these variables should be taken out of
            % the loop.
            integ_norm_sum_avg_q(ind_slow+1,:,:,:) = data_process_line.integ_norm_sum_avg_q;
        end
%         if (isfield(data_process_line,'norm_sum_avg'))
%             data_process_norm_sum_avg(:,ind_slow+1,:) = data_process_line.norm_sum_avg;
%         end
    end
end %end going through scanlines

% SASTT
% store temporary arrays from parfor loop in return structure
if (exist('integ_norm_sum_avg','var'))
    integ.norm_sum_avg = squeeze(integ_norm_sum_avg(1,:,:));
end
if (exist('integ_norm_sum_avg_q','var'))
    integ.norm_sum_avg_q = squeeze(integ_norm_sum_avg_q(1,:,:,:));
end
if (isfield(data_process,'sym_data'))
    data_process.sym_data = data_process_sym_data;
end
if (isfield(data_process,'trans_data'))
    data_process.trans_data = data_process_trans_data;
end
if (isfield(data_process,'sector_data_q'))
    data_process.sector_data_q = data_process_sector_data_q;
end
if (isfield(data_process,'sector_data'))
    data_process.sector_data = data_process_sector_data;
end


% Convert the phase of the fft coefficient to orientation in degrees
% flag for ignoring the sign of f2_phase
ignore_f2p_sign = 0;
orientation = phase2orientation(f2_phase,integ.phi_det,ignore_f2p_sign);
degree_orient = f2_amp ./ (f1_amp + 1e-6); % to avoid a division by zero
          
% Flipping data accounting for fast_axis and snake scan
% and compile return structure four_dat
[four_dat.f1_amp, ~]    = adjust_projection(f1_amp,    par.snake_scan, par.fast_axis_x, pos_data);
four_dat.f2_amp    = adjust_projection(f2_amp,    par.snake_scan, par.fast_axis_x, pos_data);
four_dat.f2_phase  = adjust_projection(f2_phase,  par.snake_scan, par.fast_axis_x, pos_data);
four_dat.orientation  = adjust_projection(orientation,  par.snake_scan, par.fast_axis_x, pos_data);
four_dat.degree_orient  = adjust_projection(degree_orient,  par.snake_scan, par.fast_axis_x, pos_data);
four_dat.I_dev     = adjust_projection(I_dev,     par.snake_scan, par.fast_axis_x, pos_data);
four_dat.I_cos_dev = adjust_projection(I_cos_dev, par.snake_scan, par.fast_axis_x, pos_data);
four_dat.I_pow     = adjust_projection(I_pow,     par.snake_scan, par.fast_axis_x, pos_data);
% peak fitting
four_dat.fit_peak.pos   = adjust_projection(fit_peak_pos  ,     par.snake_scan, par.fast_axis_x, pos_data);
four_dat.fit_peak.ampl  = adjust_projection(fit_peak_ampl ,     par.snake_scan, par.fast_axis_x, pos_data);
four_dat.fit_peak.width = adjust_projection(fit_peak_width,     par.snake_scan, par.fast_axis_x, pos_data);
% simple peak analysis
four_dat.simple_peak.pos   = adjust_projection(simple_peak_pos  ,     par.snake_scan, par.fast_axis_x, pos_data);
four_dat.simple_peak.ampl  = adjust_projection(simple_peak_ampl ,     par.snake_scan, par.fast_axis_x, pos_data);
four_dat.simple_peak.width = adjust_projection(simple_peak_width,     par.snake_scan, par.fast_axis_x, pos_data);

if (~isempty(data_adjusted))
    four_dat.positions_out = data_adjusted.positions_out;
    four_dat.scan_num      = data_adjusted.scan_num;
    four_dat.scan_point    = data_adjusted.scan_point;
end

if exist('pos_data','var')
    four_dat.original_pos_data = pos_data;
end

if (~isempty(fnames.basedir_trans_data))
    % correct for transmission data: only f1 and f2
    for ind = 1:1:size(four_dat.f1_amp,3)
        four_dat.f1_amp(:,:,ind) = four_dat.f1_amp(:,:,ind)./trans_data;
        four_dat.f2_amp(:,:,ind) = four_dat.f2_amp(:,:,ind)./trans_data;
    end
    four_dat.f1_amp(~isfinite(four_dat.f1_amp)) = 0; % Remove NaN and Inf that come from transmission = 0
    four_dat.f2_amp(~isfinite(four_dat.f2_amp)) = 0;
end

% averaged background
four_dat.I_bgr_all = I_bgr_all;

% save parameters in return structure to avoid that current settings are
% overwritten upon reloading the result of previous processing
four_dat.par = par;
four_dat.fmt = fmt;
four_dat.fnames = fnames;
four_dat.integ = integ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f1_amp_line, f2_amp_line, f2_phase_line, I_cos_dev_line, I_dev_line, I_pow_line, ...
    fit_peak_line, ...
    simple_peak_line, ...
    data_process_line] = ...
    process_one_line(par, fnames, mesh_scan_mode, file_ind, format_05, ...
    ind_slow, ind_max_slow, ind_max_fast, ...
    I_bgr_all,integ, ...
    original_trans_data)


f1_amp_line = zeros(ind_max_fast+1,length(par.r_sum));
f2_amp_line = f1_amp_line;
f2_phase_line = f1_amp_line;
I_cos_dev_line = f1_amp_line;
I_dev_line = f1_amp_line;
I_pow_line = f1_amp_line;
fit_peak_line.ampl = f1_amp_line;
fit_peak_line.pos = f1_amp_line;
fit_peak_line.width = f1_amp_line;
simple_peak_line = fit_peak_line;

% SAS TT
data_process_line = {};
if (par.prepare_segmentation)
    data_process_line.sym_data = zeros(ind_max_fast+1,size(integ.norm_sum,1));
    if (~isempty(fnames.basedir_trans_data))
        data_process_line.trans_data = zeros(ind_max_fast+1,1);
    end
end
if (par.prepare_SASTT)
    num_segments = size(integ.norm_sum, 2);
    data_process_line.sector_data = zeros(ind_max_fast+1,num_segments/2);
    data_process_line.integ_norm_sum_avg = zeros(num_segments/2,length(par.r_sum{par.which_qrange}));
    if (isfield(par,'qresolved_q') && (~isempty(par.qresolved_q)))
        data_process_line.sector_data_q = zeros(ind_max_fast+1,num_segments/2,numel(par.qresolved_q));
        data_process_line.integ_norm_sum_avg_q = zeros(num_segments/2,max(cellfun(@numel,par.qresolved_q)),numel(par.qresolved_q));
    end
end

if ((~isfield(fnames,'scan_no_step')) || (isempty(fnames.scan_no_step)))
    scan_no_step = 1;
else
    scan_no_step = fnames.scan_no_step;
end

if ~isfield(fnames, 'fext')
    fnames.fext = '.h5';
end

% online mode
waited_for_data = 0;

if (~mesh_scan_mode)
    [ filename, format_05 ] = ...
        filename_integrated_intensities(fnames.basedir_integrated_data, ...
        fnames.basename_integrated_data, ...
        file_ind,format_05,fnames.fext);
    [ filename_next, format_05 ] = ...
        filename_integrated_intensities(fnames.basedir_integrated_data, ...
        fnames.basename_integrated_data, ...
        file_ind+scan_no_step,format_05,fnames.fext);
    
    % wait for the next file to be available
    sleep_time_sec = 60;
    tic;
    if (ind_slow+scan_no_step <= ind_max_slow)
        % check that the next file is there to ensure that the current
        % one is written
        fid = utils.fopen_until_exists(filename_next,'RetryReadSleep',sleep_time_sec, ...
            'RetrySleepWhenFound',sleep_time_sec);
        fclose(fid);
        % if waiting for data was necessary then this script is running
        % online in parallel to data integration
        if (toc >= sleep_time_sec)
            waited_for_data = 1;
        end
    end
    
    % check that the current file exists
    fid = utils.fopen_until_exists(filename,'RetryReadSleep',sleep_time_sec, ...
        'RetrySleepWhenFound',sleep_time_sec);
    % in online mode wait before reading the last data file to
    % ensure that it is completely written
    if ((ind_slow == ind_max_slow) && (waited_for_data))
        pause(sleep_time_sec);
    end
    fclose(fid);

    % load radially integrated data
    if strcmpi(fnames.fext,'.h5')
        one_line = io.HDF.hdf5_load(filename);
    else
        one_line = load(filename,'I_all');
    end
end

% do it over every point in the scanline
for (ind_fast=0:ind_max_fast)
    if (ind_fast+1 > size(one_line.I_all,3))
        fprintf('%d point(s) missing\n',ind_max_fast+1-size(one_line.I_all,3));
        break;
    end
    if (mesh_scan_mode)
        I = one_line.I_all(:,:,ind_fast+1+ind_slow*(ind_max_fast+1));
    else
        I = one_line.I_all(:,:,ind_fast+1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % subtract the background
    if (~isempty(I_bgr_all))
        if (size(I_bgr_all,1) > size(I,1))
            I_bgr_all = I_bgr_all(1:size(I,1),:);
        end

        % only subtract available intensities, not the with
        % -1 flagged invalid ones
        ind_pos = intersect(find(I >= 0), find(I_bgr_all >= 0));
        I(ind_pos) = I(ind_pos) - I_bgr_all(ind_pos);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (par.prepare_segmentation)
        I_sym = sum(I.*integ.norm_sum,2)./sum(integ.norm_sum,2);
        I_sym(isnan(I_sym)) = 0;
        I_sym(I_sym < 0) = 0;
        if (isempty(fnames.basedir_trans_data))
            % in case no transmission is taken into account
            data_process_line.sym_data(ind_fast+1,:) = single(I_sym);
        else
            % correct segmentation data after transmission
            data_process_line.trans_data(ind_fast+1) = single(original_trans_data(ind_fast+1,ind_slow+1));
            data_process_line.sym_data(ind_fast+1,:) = single(I_sym./(original_trans_data(ind_fast+1,ind_slow+1)+eps));
        end
    end

    % average over the number of segments
    num_segments = size(I,2);
    if ( num_segments > 1)

        % Split the loaded intensities up in two semi-circles.
        % It will be assumed that these two semi circles are equivalent due
        % to inversion symmetry
        I1 = I(:,1:num_segments/2)';
        I9 = I(:,num_segments/2+1:num_segments)';  % Name I9 is historical from the time where 16 segments were used, hence there were 8 elements in this process
        norm_sum1 = integ.norm_sum(:,1:num_segments/2)';
        norm_sum9 = integ.norm_sum(:,num_segments/2+1:num_segments)';

        % get the deviation between valid equivalent intensities
        ind1_valid = find( I1 >= 0 );
        ind9_valid = find( I9 >= 0 );
        ind_non_zero = find(I1 + I9 ~= 0);
        ind_valid = intersect( intersect(ind1_valid,ind9_valid), ind_non_zero);
        I_diff = zeros(size(I1));
        I_diff(ind_valid) = ...
            (I1(ind_valid) - I9(ind_valid)) ./ (I1(ind_valid) + I9(ind_valid));
        I_diff = mean(I_diff,1);


        % negative intensities signal that no data are available for
        % this segment (or that too much background has been
        % subtracted) and have to be treated separately
        ind1 = find( I1 < 0 );
        ind9 = find( I9 < 0 );

        % intensities that are negative in set1 and not in set2 are
        % taken from 2
        ind = setdiff(ind1,ind9);
        I1(ind) = I9(ind);
        norm_sum1(ind) = norm_sum9(ind);
        %             figure(6); I_plot=I1; I_plot(find(I_plot<0))=1e3; imagesc(I_plot); axis xy;
        % intensities that are negative in set2 and not in set1 are
        % taken from 1
        ind = setdiff(ind9,ind1);
        I9(ind) = I1(ind);
        norm_sum9(ind) = norm_sum1(ind);

        % average the equivalent semi circles
        %             I_avg = 0.5 * (I1+I9);
        I_avg = (I1.*norm_sum1+I9.*norm_sum9)./(norm_sum1+norm_sum9+eps);
        norm_sum_avg = norm_sum1+norm_sum9;  % Number of points used in I_avg
        I_avg(norm_sum_avg==0) = -1;    % Set to -1 where no data was measured
        %             figure(6); I_plot=I_avg; I_plot(find(I_plot<0))=1e3;imagesc(I_plot); axis xy;

        % get the indices to intensities that are negative in both sets
        ind_org = intersect(ind1,ind9);
        % remove indices that have not two positive neighbouring
        % intensity values
        ind = setdiff(ind_org, ind_org-1);
        ind = setdiff(ind    , ind_org+1);
        ind_max = size(I_avg,1) * size(I_avg,2);
        ind = setdiff(ind, [1 ind_max]);
        % replace the negative intensities with the average of the
        % neighbouring pixels
        %             I_avg(ind) = 0.5 * (I_avg(ind-1)+I_avg(ind+1));
        I_avg(ind)  = ( I_avg(ind-1).*norm_sum_avg(ind-1)+I_avg(ind+1).*norm_sum_avg(ind+1) )./(norm_sum_avg(ind-1)+norm_sum_avg(ind+1)+eps);
        norm_sum_avg(ind) = norm_sum_avg(ind-1)+ norm_sum_avg(ind+1);
        %             figure(6); I_plot=I_avg;I_plot(find(I_plot<0))=1e3;imagesc(I_plot); axis xy;

        % get the indices of the so far not treated intensities
        ind_org = setdiff(ind_org,ind);
        % if there is positive intensity at one side use it
        ind = setdiff(ind_org-1, ind_org);
        ind = setdiff(ind, 0);
        I_avg(ind+1) = I_avg(ind);
        norm_sum_avg(ind+1) = norm_sum_avg(ind);
        %             figure(6); I_plot=I_avg;I_plot(find(I_plot<0))=1e3;imagesc(I_plot); axis xy;
        % get the indices of the so far not treated intensities
        ind_org = setdiff(ind_org,ind+1);
        % if there is positive intensity at the other side use it
        ind = setdiff(ind_org+1, ind_org);
        ind = setdiff(ind, ind_max+1 );
        I_avg(ind-1) = I_avg(ind);
        norm_sum_avg(ind-1) = norm_sum_avg(ind);
        %             figure(6); I_plot=I_avg;I_plot(find(I_plot<0))=1e3;imagesc(I_plot); axis xy;

        % set the remaining intensities to 0
        ind_org = setdiff(ind_org,ind-1);
        I_avg(ind_org) = 0;
        %             figure(6); I_plot=I_avg;I_plot(find(I_plot<0))=1e3;imagesc(I_plot); axis xy;

        %         % ease handling by ensuring positive values
        %         I_avg(I_avg <= 0) = 1e-19;

    else
        % One angular segment
        I_avg = I';
        norm_sum_avg = norm_sum';
        I_diff = zeros(size(I));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate matrix for q-resolved
    if isfield(par,'qresolved_q')&&(~isempty(par.qresolved_q))
        for ind = 1:numel(par.qresolved_q)
            r_curr = par.qresolved_q{ind};
            %                I_avg_q = sum(I_avg(:,r_curr), 2); % sum over the segments
            % weighted average over the segments
            % Here we can choose either 1) the average intensity per
            % pixel, 2) the total number of photons in the segment. I
            % chose (1) because it seems with value and shape closest to whatever
            % we were using before. This change probably means there
            % should be a similar change made in the SASTT
            % functions.
            I_avg_q = sum( I_avg(:,r_curr).*norm_sum_avg(:,r_curr), 2 )./( sum(norm_sum_avg(:,r_curr), 2)+eps );
            if (par.prepare_SASTT)
                data_process_line.sector_data_q(ind_fast+1,:,ind) = single(I_avg_q);
                data_process_line.integ_norm_sum_avg_q(:,1:length(r_curr), ind) = single(norm_sum_avg(:,r_curr));
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate the Fourier components
    for (ind = 1:length(par.r_sum))
        r_curr = par.r_sum{ind};

        % initialize return values
        % remarks:
        % all elemnts zero and NaN are not very aesthetic as not-performing
        % peak analysis results in an array of zeros or NaNs. 
        % But the alternative of deciding at a higher level if peak fitting
        % will be performed and carrying this through the whole code seems
        % to be even less aesthetic. 
        % The alternative of assigning an empty value [] seems to lead to
        % array element deletion and thus resizing, i.e. is highly
        % undesirable in terms of memory management and handling.
        fit_peak1.ampl = 0;
        fit_peak1.pos = 0;
        fit_peak1.width = 0;
        simple_peak1 = fit_peak1;
        
        if (((size(par.r_bgr,1) == 1) && (par.r_bgr(ind) > 0)) || ...
                ((size(par.r_bgr,1) >  1) && (par.r_bgr(ind,1) > 0)))
            % perform simple background interpolation assuming the
            % pre-defined power-of-q trend
            [ I_0, I_1, I_bgr, I_bgr_pow, r_curr_incl_bgr, r_curr_start_ind ] = interpolate_bgr(I_avg,r_curr,par,ind);

            % debugging plot to check background subtraction
            debug_plot = 0;
            if ((debug_plot) && (ind_slow == 0))
                figure(202);
                hold off;
                clf;
                r_plot_from = r_curr(1) - 50;
                if (r_plot_from < 1)
                    r_plot_from = 1;
                end
                r_plot_to = r_curr(end) + 50;
                if (r_plot_to > size(I_avg,2))
                    r_plot_to = size(I_avg,2);
                end
                r_plot = r_plot_from:r_plot_to;
                semilogy( r_plot,mean(I_avg(:,r_plot),1) );
                hold all;
                semilogy( r_curr_incl_bgr,mean(I_bgr,1) );
                title( sprintf('file %d: radius %d-%d',file_ind,r_curr(1),r_curr(end)) );
                drawnow;
            end

            if (par.peakfit(ind))
                % attempt fitting a Gaussian peak
                % Please note: If peak fitting rather than just simple peak
                % analysis is performed, then the background is often
                % estimated better than with the simple estimation
                % beforehand. Nevertheless, this better estimate is not
                % used currently, as the peak fitting operates on the
                % average of all segments, i.e. azimuthally averaged data
                % and thus the background estimation is azimuthally
                % averaged as well. This could be fixed by performing the
                % peak-fitting for each segment. This is technically
                % straight forward but would slow down the analysis and
                % poses the question if one really wants to plot the
                % multiple results. On the pro side, a peak-fit for each
                % azimuth offers new analysis options and promises better
                % results in case of peaks that are azimuthally well
                % defined, i.e. washed out by averaging over the azimuth.
                [ fit_peak1, simple_peak1, I_bgr_fit ] = fit_peak(I_avg,r_curr,r_curr_incl_bgr,I_0,I_1,I_bgr,I_bgr_pow,par,ind);
            end
        else
            % no background subtraction, i.e. zero background
            I_bgr = zeros(size(I_avg,1),length(r_curr));
            r_curr_start_ind = 1;
        end

        % Obtain intensity on average over r, i.e. q.
        % Use number of pixels that went into the calculation of each
        % individual intensity value (norm_sum) as weighing factor for
        % the calculation of this average. 
        I_avg_r = sum((I_avg(:,r_curr)-I_bgr(:,r_curr_start_ind:(r_curr_start_ind+length(r_curr)-1))) ...
            .* norm_sum_avg(:,r_curr), 2)./(sum(norm_sum_avg(:,r_curr), 2)+ eps);



        % calculate the fourier coefficients
        [f1_amp1, f2_amp1, f2_phase1, I_cos_dev1] = fourier_coefficients(I_avg_r);

        % calculate the cosine difference
        % asymmetry of the intensities
        % compares one segment to the equivalent one
        I_dev1 = mean( I_diff(r_curr) );

        % average over all segments with positive intensities, i.e.,
        % skip unavailable intensities
        %
        % Warning: I am not sure if this part should also include
        % normalization by norm_sum. Currently it does not.
        I_seg_sum = zeros(1,length(r_curr));
        for (ind2_r=1:length(r_curr))
            no_of_el = 0;
            ind2 = r_curr(ind2_r);
            for (ind1=1:size(I_avg,1))
                if (I_avg(ind1,ind2) >= 0)
                    I_seg_sum(1,ind2_r) = I_seg_sum(1,ind2_r) + I_avg(ind1,ind2);
                    no_of_el = no_of_el +1;
                end
            end
            if (no_of_el > 1)
                I_seg_sum(1,ind2_r) = I_seg_sum(1,ind2_r) / no_of_el;
            end
        end
        % use linear regression to determine the intensity power law,
        % not taking error bars into account
        ind_pos = find(I_seg_sum > 0);
        p = 0;
        if (length(ind_pos) > 1)
            p = polyfit(log(r_curr(ind_pos)),log(I_seg_sum(ind_pos)),1);
        end
        I_pow1 = p(1);

        %             if (ind > 1)
        %                 figure(1);hold off;loglog(r_curr,I_seg_sum);hold all; loglog(r_curr, r_curr.^p(1) * exp(p(2)) );
        %                 pause(1);
        %             end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %save the averaged data in the chosen q-range
        if (par.prepare_SASTT) && (ind == par.which_qrange)
            data_process_line.sector_data(ind_fast+1,:) = single(I_avg_r);
            % This does not depend on the outer-loop variable ind_fast,
            % i.e. is set more often than needed.
            data_process_line.integ_norm_sum_avg = single(norm_sum_avg(:,r_curr));
        end
        %             if (par.prepare_SASTT || par.prepare_segmentation)
        %                 if ind == par.which_qrange
        %                     counter = counter + 1;
        %                 end
        %             end

        % save return values in the corresponding arrays
        % Fourier analysis
        indy = ind_fast+1;
        f1_amp_line(indy,ind) = f1_amp1;
        f2_amp_line(indy,ind) = f2_amp1;
        f2_phase_line(indy,ind) = f2_phase1;
        % debug value
        I_cos_dev_line(indy,ind) = I_cos_dev1;
        % exp;onent of intensity decay
        I_pow_line(indy,ind) = I_pow1;
        % peak analysis / fitting
        fit_peak_line.pos(indy,ind) = fit_peak1.pos;
        fit_peak_line.ampl(indy,ind) = fit_peak1.ampl;
        fit_peak_line.width(indy,ind) = fit_peak1.width;
        simple_peak_line.pos(indy,ind) = simple_peak1.pos;
        simple_peak_line.ampl(indy,ind) = simple_peak1.ampl;
        simple_peak_line.width(indy,ind) = simple_peak1.width;
        I_dev_line(indy,ind) = I_dev1;
    end
    %         counter = counter + 1; %% Increase counter for the point number in the full 2D scan
end %end going through points in scanline
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I_all] = calc_I_point_avg(I_all_prev,fnames,plot_title)
% Warning, this function does not use norm_sum. Should be improved to
% provide an appropriate background normalization and tested.


% average over all pixels with positive intensities, i.e., skip
% negative intensities
no_of_radii = size(I_all_prev,1);
no_of_segments = size(I_all_prev,2);
% points to average, [] for all
point_range = fnames.bgr_point_range;
if (isempty(point_range))
    point_range = 1:size(I_all_prev,3);
    point_range_unfil = point_range;
end

% debug plot with the median of all data points and the selected ones
% and a second plot with the resulting background intensity profile without
% and with median filtering (if activated)
debug_plot = 0;

% if requested, return the average of the typical rather than all intensity
% lines -- selected via a very simple median-filtering
if (isfield(fnames,'bgr_max_dev_from_median') && (~isempty(fnames.bgr_max_dev_from_median)))
    if (isfield(fnames,'bgr_median_radii_range') && (~isempty(fnames.bgr_median_radii_range)))
        bgr_median_radii_range = fnames.bgr_median_radii_range;
    else
        bgr_median_radii_range = 1:no_of_radii;
    end
    % determine median intensity over specified radius-range and all
    % segments at each scan-point 
    fil_ind = 1;
    % repeat filtering until all points agree with the median, or the max.
    % no. of iterations has been reached
    while (fil_ind <= 5) 
        % filter with twice the threshold in the first iteration as the
        % median may be off due to many non-background points
        if (fil_ind == 1)
            bgr_max_dev_from_median = 2.0 * fnames.bgr_max_dev_from_median;
        else
            bgr_max_dev_from_median = fnames.bgr_max_dev_from_median;
        end
        no_of_points = numel(point_range);
        I_all = zeros(1,no_of_points);
        for (ind3=1:no_of_points)
            I_consider = I_all_prev(bgr_median_radii_range,:,point_range(ind3));
            I_all(ind3) = median(I_consider(I_consider >= 0),1);
        end
        close_to_median = abs(I_all/median(I_all) -1) <= bgr_max_dev_from_median;
        point_range = point_range(close_to_median);
        fprintf('Iteration %d: %d out of %d background points filtered out due to deviation of more than %.2f from the median intensity\n', ...
            fil_ind, no_of_points-numel(point_range),no_of_points, bgr_max_dev_from_median);
        if (fil_ind == 1)
            I_all_1 = I_all;
        end
        if (no_of_points-numel(point_range) == 0) && (fil_ind > 1)
            break;
        end
        fil_ind = fil_ind +1;
    end
    % debug plot with the median of all data points and the selected ones
    debug_plot = 1;
    if (debug_plot)
        figure(3); 
        hold off;
        clf;
        subplot(2,1,1);
        plot(point_range_unfil,I_all_1);
        hold all;
        plot(point_range,I_all(close_to_median),'x');
        legend('all points','median filtered');
        xlabel('point no.');
        ylabel('median intensity over specified range and all segments');
        title(plot_title);
    end    
end

% average intensities over scan-points
I_all = zeros(no_of_radii,no_of_segments);
for (ind1=1:no_of_radii)
    for (ind2=1:no_of_segments)
        I_all(ind1,ind2) = mean(I_all_prev(ind1,ind2,point_range(I_all_prev(ind1,ind2,point_range) >= 0)));
    end
end
% flag unavailable average due to unavailable data points as -1
I_all(isnan(I_all)) = -1;

% debug plot with filtered and unfiltered average intensity
if (debug_plot)
    % repeat unfiltered, for comparison: 
    % average intensities over scan-points
    I_all_unfil = zeros(no_of_radii,no_of_segments);
    for (ind1=1:no_of_radii)
        for (ind2=1:no_of_segments)
            I_all_unfil(ind1,ind2) = mean(I_all_prev(ind1,ind2,point_range_unfil(I_all_prev(ind1,ind2,point_range_unfil) >= 0)));
        end
    end
    % flag unavailable average due to unavailable data points as -1
    I_all_unfil(isnan(I_all)) = -1;

    figure(3); 
    subplot(2,1,2);
    hold off;
    x_data = 1:size(I_all,1);
    for (ind2 = 1:no_of_segments)
        I_plot = I_all_unfil(:,ind2);
        plot_ind = I_plot >= 0;
        semilogy(x_data(plot_ind),I_plot(plot_ind),'b');
        if (ind2 == 1)
            hold all;
        end
        I_plot = I_all(:,ind2);
        plot_ind = I_plot >= 0;
        semilogy(x_data(plot_ind),I_plot(plot_ind),'g');
    end
    xlabel('radius [pixel]');
    ylabel('average intensity [a.u.]');
    legend('unfiltered','filtered');
    title(sprintf('Intensities for all %d segments',no_of_segments));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ filename, format_05 ] = ...
    filename_integrated_intensities(basedir_integrated_data, ...
    basename_integrated_data, ...
    scan_no,format_05, fext)


if (~isempty(basename_integrated_data))
    if (contains(basename_integrated_data,'%'))
        % the format has been specified
        filename = [ basedir_integrated_data ...
            sprintf(basename_integrated_data,scan_no) ];
    else
        % try to get the format right
        filename = [ basedir_integrated_data basename_integrated_data ...
            num2str(scan_no,'%05d') '_00000_00000_integ' fext ];
        if (~exist(filename,'file'))
            filename = [ basedir_integrated_data basename_integrated_data ...
                num2str(scan_no,'%05d') '_00000_integ' fext ];
            if (~exist(filename,'file'))
                fprintf('%s not found. Trying normal scan convention.\n',filename);
                filename = [ basedir_integrated_data basename_integrated_data ...
                    num2str(scan_no,'%05d') '_00000_00_integ' fext ];
                if (~exist(filename,'file'))
                    error([ filename ' not found']);
                end
            end
        end
    end
else
    % older file name formats
    if ((format_05 == -1) || (format_05 == 1))
        filename = [ basedir_integrated_data 'image_1_' ...
            num2str(scan_no,'%05d') ...
            '_00000_00_00000_integ' fext ];
        if (format_05 == -1)
            if (exist(filename,'file'))
                format_05 = 1;
            else
                format_05 = 0;
            end
        end
    end
    if (format_05 == 0)
        filename = [ basedir_integrated_data 'image_1_' ...
            num2str(scan_no,'%04d') ...
            '_00000_00_00000_integ' fext ];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ I_0, I_1, I_bgr, I_bgr_pow, r_curr_incl_bgr, r_curr_start_ind ] = interpolate_bgr(I_avg,r_curr,par,ind)


% if only one number is specified then the interpolation
% range for the background startsin direct proximity to the
% summation range for the intensity
%
% Warning, this interpolation does not take into account
% the number of elements used in the calcualtion of average
% intensities (norm_sum) I cannot include this
% normalization without testing with extraction of an
% actual peak because norm_sum changes would modify the
% trend of the data vs q.
if ((ind <= size(par.r_bgr,1)) && (2 == size(par.r_bgr,2)))
    r_bgr_far = par.r_bgr(ind,1);
    r_bgr_close = par.r_bgr(ind,2);
else
    % default is to start at the neighboring pixel
    r_bgr_far =   par.r_bgr(ind);
    r_bgr_close = 1;
end

% get indices and center value of the background patches on
% both sides of the current region
r_bgr_left  = (r_curr(1)-r_bgr_far):(r_curr(1)-r_bgr_close);
r_left = 0.5 * (r_bgr_left(1) + r_bgr_left(end));
r_bgr_right = (r_curr(end)+r_bgr_close):(r_curr(end)+r_bgr_far);
r_right = 0.5 * (r_bgr_right(1) + r_bgr_right(end));
% index over the background including peak range
r_curr_incl_bgr = r_bgr_left(1):r_bgr_right(end);
% start index of r_curr within r_curr_incl_bgr
r_curr_start_ind = r_bgr_far +1;
% calculate the average intensity on the left and right
% side of the current region
I_bgr_left = sum( I_avg(:,r_bgr_left), 2) / length(r_bgr_left);
I_bgr_right = sum( I_avg(:,r_bgr_right), 2) / length(r_bgr_right);
% user specified intensity decay as background
I_bgr_pow = par.I_bgr_pow(ind);
if (I_bgr_pow >= 0)
    % interpolate a I0 + I1 * q^-I_bgr_pow background
    I_pow_div = ( r_left^(-I_bgr_pow) - r_right^(-I_bgr_pow) );
    if (abs(I_pow_div) > 1e-30)
        I_1 = (I_bgr_left - I_bgr_right) / I_pow_div;
        I_0 = I_bgr_right - I_1 * r_right^(-I_bgr_pow);
    else
        % in cases like I_bgr_pow == 0
        I_1 = zeros(size(I_bgr_left));
        I_0 = (I_bgr_left + I_bgr_right) ./ 2;
    end
    % return as well background over the current region
    I_bgr = I_0 + I_1 * r_curr_incl_bgr.^(-I_bgr_pow);
else
    % alternatively assume that the intensity decay power law leads to
    % background on left and right side being equal, i.e. determine
    % exponent of intensity decay based on this assumption
    I_bgr_pow = -log(I_bgr_left./I_bgr_right) / log(r_left/r_right);
    % check that the intensity power law exponent is finite and 
    % constrain to the range expected: [0 4] (or a wider range)
    I_bgr_pow(isinf(I_bgr_pow)) = 0;
    I_bgr_pow(isnan(I_bgr_pow)) = 0;
    I_bgr_pow(I_bgr_pow < -0.1) = -0.1;
    I_bgr_pow(I_bgr_pow >  4.1) =  4.1;
    I_1 = (I_bgr_left - I_bgr_right) ./ (r_left.^(-I_bgr_pow) - r_right.^(-I_bgr_pow));
    I_1(isinf(I_1)) = 0;
    I_0 = I_bgr_left - I_1 .* r_left.^(-I_bgr_pow);
    % return as well background over the current region
    I_bgr = I_0 + I_1 .* r_curr_incl_bgr.^(-I_bgr_pow);
end

% debugging code
debug_plot = 0;
if (debug_plot)
    figure(201);
    hold off;
    plot(r_curr_incl_bgr,mean(I_avg(:,r_curr_incl_bgr),1));
    hold on;
    plot(r_curr_incl_bgr,mean(I_bgr,1),'.');
    legend('intensity','background');
    title(sprintf('I\\_bgr\\_pow = %.1f',mean(I_bgr_pow)));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ fit_peak_one, simple_peak_one, I_bgr ] = ...
    fit_peak(I_avg,r_curr,r_curr_incl_bgr,I_0,I_1,I_bgr,I_bgr_pow,par,ind)

% fit to the average of all azimuthal segments
% and use the background estimation as starting guess for the fit
I_for_fit = mean(I_avg(:,r_curr_incl_bgr));
I0 = mean(I_0);
I1 = mean(I_1);
Ibgr = mean(I_bgr,1);

% simple peak interpretation based on intensity above simple background
% interpolation
I_minus_bgr = I_for_fit - Ibgr;

% initialize return values
simple_peak_one.ampl = 0;
simple_peak_one.pos = 0;
simple_peak_one.width = 0;
fit_peak_one.pos = 0;
fit_peak_one.ampl = 0;
fit_peak_one.width = 0;

% maximum intensity position as starting guess for peak position
[ peak_int, peak_pos ] = max(I_minus_bgr);
peak_pos = peak_pos(1);

% no intensity: return initialization values
if (isnan(peak_int))
    return
end

% go to lower index side until half intensity reached
curr_pos = peak_pos;
while ((I_minus_bgr(curr_pos) > 0.5*peak_int) && (curr_pos > 1))
    curr_pos = curr_pos - 1;
end
% interpolate half-intensity index on this side
if (curr_pos < length(I_minus_bgr))
    if (I_minus_bgr(curr_pos+1) == I_minus_bgr(curr_pos))
        left_pos = curr_pos;
    else
        left_pos = curr_pos + (0.5*peak_int - I_minus_bgr(curr_pos)) / (I_minus_bgr(curr_pos+1) - I_minus_bgr(curr_pos));
    end
else
    left_pos = curr_pos -1;
end
% go to upper index side until half intensity reached
curr_pos = peak_pos;
while ((I_minus_bgr(curr_pos) > 0.5*peak_int) && (curr_pos < length(I_minus_bgr)))
    curr_pos = curr_pos + 1;
end
% interpolate half-intensity index on this side
if (curr_pos > 1)
    if (I_minus_bgr(curr_pos) == I_minus_bgr(curr_pos-1))
        right_pos = curr_pos;
    else
        right_pos = curr_pos + (0.5*peak_int - I_minus_bgr(curr_pos)) / (I_minus_bgr(curr_pos) - I_minus_bgr(curr_pos-1));
    end
else
    right_pos = length(I_minus_bgr);
end

% intensities, mainly for debug plot
if (left_pos < 1)
    left_pos = 1;
    I_left = I_for_fit(1);
    I_left_no_bgr = I_minus_bgr(1);
else
    if (left_pos >= length(I_for_fit))
        left_pos = length(I_for_fit);
        I_left = I_for_fit(end);
        I_left_no_bgr = I_minus_bgr(end);
    else
        I_left = I_for_fit(floor(left_pos)) + ...
            (I_for_fit(floor(left_pos) +1) - I_for_fit(floor(left_pos))) * ...
            (left_pos - floor(left_pos));
        I_left_no_bgr = I_minus_bgr(floor(left_pos)) + ...
            (I_minus_bgr(floor(left_pos) +1) - I_minus_bgr(floor(left_pos))) * ...
            (left_pos - floor(left_pos));
    end
end
if ((right_pos+1 >= length(I_for_fit)) || isnan(right_pos))
    right_pos = length(I_for_fit);
    I_right = I_for_fit(end);
    I_right_no_bgr = I_minus_bgr(end);
else
    if (right_pos <= 1)
        right_pos = 1;
        I_right = I_for_fit(1);
        I_right_no_bgr = I_minus_bgr(1);
    else
        I_right = I_for_fit(floor(right_pos)) + ...
            (I_for_fit(floor(right_pos) +1) - I_for_fit(floor(right_pos))) * ...
            (right_pos - floor(right_pos));
        I_right_no_bgr = I_minus_bgr(floor(right_pos)) + ...
            (I_minus_bgr(floor(right_pos) +1) - I_minus_bgr(floor(right_pos))) * ...
            (right_pos - floor(right_pos));
    end
end
peak_pos = 0.5 * (left_pos + right_pos);
if (peak_pos < 1)
    I_bgr_at_peak = Ibgr(1);
    I_peak = I_for_fit(1);
else
    if (peak_pos+1 > length(Ibgr))
        I_bgr_at_peak = Ibgr(end);
        I_peak = I_for_fit(end);
    else
        % interpolated background intensity at peak position
        I_bgr_at_peak = Ibgr(floor(peak_pos)) + ...
            (Ibgr(floor(peak_pos) +1) - Ibgr(floor(peak_pos))) * ...
            (peak_pos - floor(peak_pos));
        % mean background intensity in the central FWHM of the peak, use
        % with care as the background is non-linear
        I_mean_bgr_peak_center = mean(Ibgr(floor(left_pos):ceil(right_pos)));
        % interpolated peak intensity over background, can be used for caluclating the return value as well
        I_peak = I_for_fit(floor(peak_pos)) + ...
            (I_for_fit(floor(peak_pos) +1) - I_for_fit(floor(peak_pos))) * ...
            (peak_pos - floor(peak_pos));
    end
end

% convert to r_curr_incl_bgr index range
left_pos = left_pos + r_curr_incl_bgr(1) -1;
right_pos = right_pos + r_curr_incl_bgr(1) -1;

% refine peak-center estimate as average between lower and upper
% half-intensity position (and convert to r_curr_incl_bgr index range)
peak_pos = 0.5 * (left_pos + right_pos);

% compile estimates for amplitude, position and half-width of the peak
% Return average of maximum and interpolated intensity at peak-position above background as peak
% amplitude: average of peak_int and I_peak - I_bgr_at_peak. 
% Alternatively one could return the maximum intensity above background:
% peak_int.
% Or combinations thereof. 
simple_peak_one.ampl = 0.5 * (I_peak - I_bgr_at_peak + peak_int); 
simple_peak_one.pos = peak_pos;
simple_peak_one.width = 0.5 * (right_pos - left_pos);

% hard-coded debug plot
fit_debug_plot = 0;
debug_plot_shown = 0;

% debug plot, adapt thresholds to your needs
if (fit_debug_plot) && ...
        (simple_peak_one.ampl > 0.1 * I_mean_bgr_peak_center > 0.1) && ...
        (simple_peak_one.ampl > 0.1)
    debug_plot_shown = 1;
    
    if ishandle(201)
        % if the figure exists, select it
        figure(201);
    else
        % if the figure does not exist yet, create it with default size at
        % default location and double its vertical size then
        figure(201);
        figpos = get(gcf,'Position');
        figpos(2) = max([figpos(2)-figpos(4) 10]);
        figpos(4) = 2 * figpos(4);
        set(gcf,'Position',figpos);
    end
    hold off;
    clf;
    subplot(2,1,1);
    % plot data
    plot(r_curr_incl_bgr,I_for_fit,'b');
    hold on;
    % plot estimated background
    plot(r_curr_incl_bgr,Ibgr,'k');
    legend('intensity','background');
    % mark peak position
    plot([simple_peak_one.pos, simple_peak_one.pos],[0.5*(I_left+I_right), simple_peak_one.ampl + I_bgr_at_peak],'*-g');
    % mark peak width
    plot([left_pos, right_pos], [I_left, I_right],'*-g');
    title(sprintf('simple peak analysis debug plot, ind=%d',ind));
    
    subplot(2,1,2);
    % plot data after background subtraction
    plot(r_curr_incl_bgr,I_minus_bgr);
    hold on;
    % mark peak position
    plot([simple_peak_one.pos, simple_peak_one.pos],[0.5*(I_left_no_bgr+I_right_no_bgr), simple_peak_one.ampl],'*-g');
    % mark peak width
    plot([left_pos, right_pos], [I_left_no_bgr, I_right_no_bgr],'*-g');
    title(sprintf('after background subtraction, peak_{ampl}=%.3f, peak_{pos}=%.1f',simple_peak_one.ampl,simple_peak_one.pos));
end

if (par.peakfit(ind) > 1)
    % use a Gaussian peak and a power-of-q background as fit function, 
    % use the results of the simple peak analysis above as starting guess
    a_start = simple_peak_one.ampl;
    b_start = simple_peak_one.pos;
    c_start = simple_peak_one.width;
    % use the above determined background as starting guess
    % for the background
    d_start = max([0, I1]);
    % use the pre-defined or determined basckground power law
    e_start = max([ 0 mean(I_bgr_pow)]);
    f_start = max([0, I0]);
    % store starting parameters and fit function in variables
    fit_opt = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',     [                        0.0,   r_curr(1),                     1,                      0.0,       0,                    0.0  ], ...
            'Upper',     [  1.5*max([I_for_fit 1e-7]), r_curr(end), r_curr(end)-r_curr(1),  max([1.5*abs(I1) 1e-7]),       4, 1.5*max([abs(I0) 1e-7]) ], ...
            'Startpoint',[                    a_start,     b_start,               c_start,                  d_start, e_start,                 f_start ] );
    fit_type = fittype('a*exp(-((x-b)/c)^2)+d*x^(-e)+f');

    % finally optimize the parameters, i.e. perform the peak fitting
    [fit_obj,~] = fit(r_curr_incl_bgr',I_for_fit',fit_type,fit_opt);

    % compile the return values for amplitude, position and half-width of the peak
    fit_peak_one.ampl = fit_obj.a;
    fit_peak_one.pos = fit_obj.b;
    fit_peak_one.width = fit_obj.c;
    % compile return values for the background intensity, 
    % return the same background vector for each segment (!)
    I_bgr = ones(size(I_avg,1),1) * (fit_obj.d * r_curr_incl_bgr.^(-fit_obj.e) + fit_obj.f);

    % debug information
    if (fit_debug_plot)
        if (debug_plot_shown ~= 0)
            % plot resulting Gaussian
            subplot(2,1,1);
            plot(fit_obj);
        end
        fprintf('a amplitude  = %10.3f / %10.3f\nb position   = %10.3f / %10.3f\nc half-width = %10.3f / %10.3f\nd bgr. scale = %10.3e / %10.3e\ne bgr. power = %10.3f / %10.3f\nf bgr. const = %10.3f / %10.3f\n\n',...
            a_start,fit_obj.a, b_start,fit_obj.b, c_start,fit_obj.c, d_start,fit_obj.d, e_start,fit_obj.e, f_start,fit_obj.f);
    end
else
    % signal that no peak fitting has been performed
    fit_peak_one.pos = 0;
    fit_peak_one.ampl = 0;
    fit_peak_one.width = 0;
    % compile return values for the background intensity, 
    % return the same background vector for each segment (!)
    % (to be consistent with the fit to a Gaussian function)
    I_bgr = ones(size(I_avg,1),1) * (I1 * r_curr_incl_bgr.^(-mean(I_bgr_pow)) + I0);
end

if (fit_debug_plot)
    hold off;
end
