% template_fluo_tomo.m
%                                                                                                   
% to process fluorescence data from many projections
% and make tomographic reconstruction

% This code and subroutines are part of a continuous development. There is no 
% liability on PSI or cSAXS. Before publishing results using this code 
% please notify M. Guizar-Sicairos, (mguizar@gmail.com)
% or another PSI cSAXS staff member.

close all
%% Prepare substructures for paths, processing and detector calibration 
%  this section is likely to be changed only once in the beginning of the 
%  experiment.

% paths substructure
fls.    path.base_path = '../../../online_data/';
fls.    path.matlab = '../../training_ana/';
fls.    path.xray_database = './+fluo/xraydb.mat'; % this is the database with X-ray emission and absorption lines
fls.    path.output = '../../analysis_fluo/';
fls.    path.marked_lines = []; % path to file with fluorescence lines that have been marked and identified on the spectrum
fls.    path.calibration = '../../training_ana/base/+fluo/detector_calibration_200915_163627.mat'; % file with detector calibration
fls.    path.ptycho_alignment = '/das/work/p18/p18043/online_data/analysis_tomo/tomo_id_223_chip_S00328_to_S00434__recons_/stack_object_external.mat';
fls.    path.fluo_datafile = '../../analysis_fluo/tomo_id223_fluo_data_raw.mat'; % the file with decoded fluo data. 
                                                                                 % Leave empty [] if no decoded fluo data available

% processing substructure
fls.    proc.metadata_source = 'orchestra'; % options are 'spec', 'orchestra' and 'sgalil'
fls.    proc.file_type = 'raw'; % options are 'nexus' and 'raw'
fls.    proc.motorx = []; % x-axis motor NOTE: ignored if the metadata_source = 'sgalil' or 'orchestra'
fls.    proc.motory = []; % y-axis motor NOTE: ignored if the metadata_source = 'sgalil' or 'orchestra'

fls.    proc.dead_time_correction = false; % as default better always correct the dead time
fls.    proc.stitched = 'none'; % options are 'none', 'blocks', 'lines'
fls.    proc.use_interpolation = false; % option for using interpolation when mapping the fluorescence data
fls.    proc.mapping_mode = true; % if mapping mode is true, then 2D map is assumed, otherwise just multiple spectra like in loopscan

% other parameters
fls.    calib.scan_num = [327]; % scans for calibration, can be a range of scans or a single scan
                                 % range of scans is summed up in a single spectrum, but if the scans represent
                                 % different projections instead of the stitches, only one scan will be
                                 % displayed during the calibration
fls.    calib.tested = false;

addpath(fls.path.matlab)
addpath('../tomo/')
fls     = fluo.unpack_xray_database(fls);

%% Load or perform detector calibration & select peaks of interest
%  This section allows you to calibrate the detector from channel numbers
%  into the units of keV. If the fls.calib.file is provided then you will
%  be invited to inspect whether the detector is properly calibrated. 
%  However, if there is no fls.calib.file provided then you will be asked 
%  to assist the script in calibrating the spectrum summed from the scans
%  given by fls.calib.scan_num.

fls = fluo.calibrate_detector(fls);
fls = fluo.test_detector_calibration(fls);
POIs = fluo.mark_peaks_of_interest(fls);

use_auto_ID = false; % if you tested the detector calibration and found that it failed -> set use_auto_ID to false
if use_auto_ID
    disp('Performing automatic peaks identification.')
    marked_peaks = fluo.auto_line_id(fls, POIs, fls.calib.energy);
else
    disp('Loading the data for manual peak identification.')
    fluo_data = fluo.fluo_read(fls.path.base_path, 'scanNr', fls.calib.scan_num, 'fluo_structure', fls); 
    marked_peaks = fluo.manual_line_id(POIs, fls.calib.energy, fluo_data);
end

%% Inspect data
%  This section allows you to load specific scan and inspect the energy
%  lines that you just identified
scans_to_inspect = [340];
fluo_data = fluo.fluo_read(fls.path.base_path, 'scanNr', scans_to_inspect, 'fluo_structure', fls);
fluo.display_peaks(fls, marked_peaks, fluo_data)

%% Load alignment dataset
alignment = load(fls.path.ptycho_alignment);
theta = alignment.theta;
total_shift = alignment.total_shift;
scanstomo = alignment.par.scanstomo;
ptycho_pixel_size = alignment.par.pixel_size;

%% Load fluorescence data
% since loading fluorescence data is a rather slow process due to the
% decoding of the raw data, in this section we have an option of saving
% the dataset after loading it. Then you can give the path to the saved
% file in the first section of this template.
energy = fls.calib.energy;
fluo_pixel_size = 100e-9; %!!! HARD CODED VALUE NEED TO GET SOMEWHERE IN SPEC
save_file_under_name = []; % if you provide the filename then the code will export loaded data for faster loading later
                           % the name can be like that: '../../analysis_fluo/tomo_id223_fluo_data_raw.mat'

if ~isempty(fls.path.fluo_datafile)
    disp('Loading fluo .mat data file.')
    fluo_data_loaded = load(fls.path.fluo_datafile);
    fluo_data = fluo_data_loaded.fluo_data;
else
    disp('Loading raw fluorescence data (i.e. not decoded raw data).')
    fluo_data = fluo.fluo_read(fls.path.base_path, 'scanNr', scanstomo, 'fluo_structure', fls);
    if ~isempty(save_file_under_name)
        save(save_file_under_name,'fluo_data','scanstomo', 'theta', 'total_shift', 'fluo_pixel_size', '-v7.3','-nocompression')
    end
end

%% Inspect the data
% Here you can inspect the loaded projections, for this you need to
% select the number of the peak that you want to inspect (assuming that 
% you marked peaks of interest in one of the previous steps)

peak_number = 11; % to list available peaks just type marked_peaks in the command window

ind_cen = marked_peaks.index_center(peak_number);
ind_min = ind_cen - floor(marked_peaks.bandwidth(peak_number) / 2);
ind_max = ind_cen + floor(marked_peaks.bandwidth(peak_number) / 2);
ch_range = [ind_min:ind_max];

[num_ch, size_y, size_x] = size(fluo_data{1});
padded_size_x = ceil(size_x * 1.1);
padded_size_y = ceil(size_y * 1.1);
stack_object = zeros(padded_size_y, padded_size_x, numel(scanstomo));

f = waitbar(0,'Arranging the projections for inspection...');
for jj = 1:numel(scanstomo)
    waitbar(jj/numel(scanstomo), f, 'Arranging the projections for inspection...');
    intermid_data = fluo_data{jj};
    intermid_data(isnan(intermid_data))=0;
    intermid_data(isinf(intermid_data))=0;
    temp_data = squeeze(sum(intermid_data(ch_range, :, :), 1));
    stack_object(:, :, jj) = utils.crop_pad(temp_data, [padded_size_y padded_size_x]);
end
close(f);
[~, ind_sort] = sort(theta);
plotting.imagesc3D(stack_object(:, :, ind_sort)), axis image xy, colorbar, colormap bone, colorbar

%% Inspect the data after applying the alignment
peak_number = 11;
ind_cen = marked_peaks.index_center(peak_number);
ind_min = ind_cen - floor(marked_peaks.bandwidth(peak_number) / 2);
ind_max = ind_cen + floor(marked_peaks.bandwidth(peak_number) / 2);

ch_range = [ind_min:ind_max];

pixel_ratio = ptycho_pixel_size / fluo_pixel_size;
total_shift_x = total_shift(:, 1) * pixel_ratio;
total_shift_y = total_shift(:, 2) * pixel_ratio;
constant_offset_x = 0;
constant_offset_y = 0;

[num_ch, size_y, size_x] = size(fluo_data{1});
padded_size_x = size_x * 2;
padded_size_y = size_y * 2;
stack_object = zeros(padded_size_y, padded_size_x, numel(scanstomo));

f = waitbar(0,'Aligning projections...');
for jj = 1:numel(scanstomo)
    waitbar(jj/numel(scanstomo), f, 'Aligning projections...');
    intermid_data = fluo_data{jj};
    intermid_data(isnan(intermid_data))=0;
    intermid_data(isinf(intermid_data))=0;
    temp_data = squeeze(sum(intermid_data(ch_range, :, :), 1));
    
    temp_data = squeeze(utils.imshift_fft(temp_data, -total_shift_x(jj)-constant_offset_x, -total_shift_y(jj)-constant_offset_x, true));
    stack_object(:, :, jj) = utils.crop_pad(temp_data, [padded_size_y padded_size_x]);
end
close(f);
[~, ind_sort] = sort(theta);
plotting.imagesc3D(stack_object(:, :, ind_sort)), axis image xy, colormap bone, colorbar

%% Perform and inspect the reconstruction
peak_number = 6;
ind_cen = marked_peaks.index_center(peak_number);
ind_min = ind_cen - floor(marked_peaks.bandwidth(peak_number) / 2);
ind_max = ind_cen + floor(marked_peaks.bandwidth(peak_number) / 2);

ch_range = [ind_min:ind_max];

pixel_ratio = ptycho_pixel_size / fluo_pixel_size;
total_shift_x = total_shift(:, 1) * pixel_ratio;
total_shift_y = total_shift(:, 2) * pixel_ratio;
constant_offset_x = 0;
constant_offset_y = 0;

[num_ch, size_y, size_x] = size(fluo_data{1});
padded_size_x = size_x * 2;
padded_size_y = size_y * 2;
stack_object = zeros(padded_size_y, padded_size_x, numel(scanstomo));

f = waitbar(0,'Aligning projections and reconstructing...');
for jj = 1:numel(scanstomo)
    waitbar(jj/numel(scanstomo), f, 'Aligning projections and reconstructing...');
    intermid_data = fluo_data{jj};
    intermid_data(isnan(intermid_data))=0;
    intermid_data(isinf(intermid_data))=0;
    temp_data = squeeze(sum(intermid_data(ch_range, :, :), 1));
    
    temp_data = squeeze(utils.imshift_fft(temp_data, -total_shift_x(jj)-constant_offset_x, -total_shift_y(jj)-constant_offset_x, true));
    stack_object(:, :, jj) = utils.crop_pad(temp_data, [padded_size_y padded_size_x]);
end
close(f);
[~, ind_sort] = sort(theta);
volume_in = stack_object;
theta_use = theta(ind_sort);
volume_in = volume_in(:, :, ind_sort);

volume_out = fluo.tomo_recon_loc(volume_in, theta_use, alignment);
plotting.imagesc3D(squeeze(volume_out)), axis image xy, colormap bone, colorbar

%% Reconstruct marked energy lines and save the resulted tomograms
[~, ind_sort] = sort(theta);
theta_use = theta(ind_sort);

filename_suffix = '_run_1';
scans_string = sprintf('tomo_fluo_scans_S%05d_to_S%05d', min(scanstomo), max(scanstomo));
fluo_tomograms = [];

for kk = 1:numel(marked_peaks.index_center)
    peak_number = kk;
    ind_cen = marked_peaks.index_center(peak_number);
    ind_min = ind_cen - floor(marked_peaks.bandwidth(peak_number) / 2);
    ind_max = ind_cen + floor(marked_peaks.bandwidth(peak_number) / 2);

    ch_range = [ind_min:ind_max];

    [num_ch, size_y, size_x] = size(fluo_data{1});
    stack_object = zeros(padded_size_y, padded_size_x, numel(scanstomo));

    f = waitbar(0,'Aligning projections and reconstructing...');
    for jj = 1:numel(scanstomo)
        waitbar(jj/numel(scanstomo), f, 'Aligning projections and reconstructing...');
        intermid_data = fluo_data{jj};
        intermid_data(isnan(intermid_data))=0;
        intermid_data(isinf(intermid_data))=0;
        temp_data = squeeze(sum(intermid_data(ch_range, :, :), 1));

        temp_data = squeeze(utils.imshift_fft(temp_data, -total_shift_x(jj)-constant_offset_x, -total_shift_y(jj)-constant_offset_x, true));
        stack_object(:, :, jj) = utils.crop_pad(temp_data, [padded_size_y padded_size_x]);
    end
    close(f);
    
    volume_in = stack_object;
    volume_in = volume_in(:, :, ind_sort);
    volume_out = fluo.tomo_recon_loc(volume_in, theta_use, alignment);
    fluo_tomograms.tomo{kk} = volume_out;
    fluo_tomograms.line_id{kk} = marked_peaks.line_id{kk};
end
filename = fullfile([fls.path.output '/' scans_string filename_suffix '.mat']);
save(filename,'fluo_tomograms', '-v7.3','-nocompression')

%% Align whole fluorescence dataset (i.e. all channels) and create a hyper stack object
% note that this is a very slow process but it is required for
% reconstructing the tomogram channel by channel
% also note that the projections of this hyper stack will be sorted by angle
[num_ch, size_y, size_x] = size(fluo_data{1});
hyper_stack_object = zeros(padded_size_y, padded_size_x, num_ch, numel(scanstomo));

f = waitbar(0,'Aligning projections for all channels. This is a slow operation...');
s = clock;
for jj = 1:numel(scanstomo)
    intermid_data = fluo_data{ind_sort(jj)};
    intermid_data(isnan(intermid_data))=0;
    intermid_data(isinf(intermid_data))=0;
    for kk = 1:num_ch
        temp_data = squeeze(intermid_data(kk, :, :));
        temp_data = squeeze(utils.imshift_fft(temp_data, -total_shift_x(ind_sort(jj))-constant_offset_x, -total_shift_y(ind_sort(jj))-constant_offset_y, true));
        hyper_stack_object(:, :, kk, jj) = utils.crop_pad(temp_data, [padded_size_y padded_size_x]);
    end
    if jj ==1
      is = etime(clock,s);
      esttime = is * numel(scanstomo);
    end
    message_out = sprintf('Aligning projections for all channels. Remaining time=%s sec', num2str(esttime-etime(clock,s), '%4.1f'));
    waitbar(jj/numel(scanstomo), f, message_out);
end
close(f)
% here we sum up the whole spectrum and reconstruct the tomogram for a
% quick inspection
stack_object = squeeze(sum(hyper_stack_object, 3));
plotting.imagesc3D(stack_object), axis image xy, colormap bone, colorbar
volume_out = fluo.tomo_recon_loc(stack_object, theta_use, alignment);
plotting.imagesc3D(volume_out), axis image xy, colormap bone, colorbar

%% Reconstruct hyper tomogram (i.e. tomogram with an additional dimension
%  represented by the channels of the detector)
filename_suffix = '_run_1';
scans_string = sprintf('hypertomo_fluo_full_scans_S%05d_to_S%05d', min(scanstomo), max(scanstomo));

volume_in =  squeeze(hyper_stack_object(:, :, 1, :));
volume_out = fluo.tomo_recon_loc(volume_in, theta_use, alignment);
size_vol = size(volume_out);

hyper_tomo = zeros(size_vol(1), size_vol(2), size_vol(3), num_ch);
hyper_tomo(:, :, :, 1) = volume_out;

f = waitbar(0,'Reconstructing tomograms for each channel. This is a slow operation...');
s = clock;
for jj = 2 : num_ch
   volume_in =  squeeze(hyper_stack_object(:, :, jj, :));
   volume_out = fluo.tomo_recon_loc(volume_in, theta_use, alignment);
   hyper_tomo(:, :, :, jj) = volume_out;
   if jj ==2
      is = etime(clock,s);
      esttime = is * num_ch;
   end
   message_out = sprintf('Reconstructing tomograms for each channel. Remaining time=%s sec', num2str(esttime-etime(clock,s), '%4.1f'));
   waitbar(jj/num_ch, f, message_out);
end
close(f);

filename = fullfile([fls.path.output '/' scans_string filename_suffix '.mat']);
save(filename,'hyper_tomo', '-v7.3','-nocompression')

%% Inspect spectrum of the hyper tomogram
spectrum_tomo = squeeze(sum(hyper_tomo, 1));
spectrum_tomo = squeeze(sum(spectrum_tomo, 1));
spectrum_tomo = squeeze(sum(spectrum_tomo, 1));
spectrum_raw = squeeze(sum(hyper_stack_object, 1));
spectrum_raw = squeeze(sum(spectrum_raw, 1));
spectrum_raw = squeeze(sum(spectrum_raw, 2));
fsp(1) = semilogy(energy, spectrum_tomo, 'b-'), hold on
fsp(2) = semilogy(energy, spectrum_raw, 'r-'), hold on
legend(fsp, {'hyper tomo spectrum', 'raw hyper stack spectrum'});
title('Summed spectrum from hyper tomo')
xlabel('Energy, [keV]')

%% Extract marked energy lines and save the resulted tomograms
[~, ind_sort] = sort(theta);
theta_use = theta(ind_sort);

filename_suffix = '_run_1';
scans_string = sprintf('hypertomo_fluo_marked_peaks_scans_S%05d_to_S%05d', min(scanstomo), max(scanstomo));
fluo_tomograms = [];

f = waitbar(0,'Saving marked lines from hyper tomo...');
for kk = 1:numel(marked_peaks.index_center)
    
    peak_number = kk;
    ind_cen = marked_peaks.index_center(peak_number);
    ind_min = ind_cen - floor(marked_peaks.bandwidth(peak_number) / 2);
    ind_max = ind_cen + floor(marked_peaks.bandwidth(peak_number) / 2);

    ch_range = [ind_min:ind_max];

    fluo_tomograms.tomo{kk} = squeeze(sum(hyper_tomo(:, :, :, ch_range), 4));
    fluo_tomograms.line_id{kk} = marked_peaks.line_id{kk};
    waitbar(kk/numel(marked_peaks.index_center), f, 'Saving marked lines from hyper tomo...');
end
close(f);
filename = fullfile([fls.path.output '/' scans_string filename_suffix '.mat']);
save(filename,'fluo_tomograms', '-v7.3','-nocompression')


plotting.imagesc3D(fluo_tomograms.tomo{3}), axis image xy, colormap bone, colorbar

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS matlab package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.