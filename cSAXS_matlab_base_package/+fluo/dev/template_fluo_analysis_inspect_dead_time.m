% template_fluo_analysis_inspect_dead_time.m  
%                                                                                                   
% template for inspecting dead time values in the fluo maps

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
fls.    path.matlab = '../../cxs_software_fluo/';
fls.    path.xray_database = './+fluo/xraydb.mat'; % this is the database with X-ray emission and absorption lines
fls.    path.probe = ''; % path to the file with reconstructed probe, only used when deconvolution is done
fls.    path.output = '../../analysis_fluo/';
fls.    path.marked_lines = []; % path to file with fluorescence lines that have been marked and identified on the spectrum

addpath(fls.path.matlab)
fls     = fluo.unpack_xray_database(fls);

% processing substructure
fls.    proc.metadata_source = 'orchestra'; % options are 'spec' and 'orchestra'
fls.    proc.file_type = 'raw'; % options are 'nexus' and 'raw'
fls.    proc.dead_time_correction = true; % as default better always correct the dead time

%% CELL
%% Load alignment dataset
alignment = load('/das/work/p17/p17648/online_data/analysis_tomo/tomo__test_id_236_S1_2_S01110_to_S01223__recons_/stack_object_external.mat');
theta = alignment.theta;
total_shift = alignment.total_shift;
scanstomo = alignment.par.scanstomo;
ptycho_pixel_size = alignment.par.pixel_size;

%% Process data
[~, ind_sort] = sort(theta);
sorted_scans = scanstomo(ind_sort);
sorted_theta = theta(ind_sort);
range = [722:775];
subfolder = 'cell_dead_time/';
if ~exist([fls.path.output subfolder], 'dir')
   mkdir([fls.path.output subfolder])
end

for jj = 1:numel(sorted_scans)
    scans_to_stich = sorted_scans(jj);

    [fluo_data_corrected, real_time_data, live_time_data, correction_coeff_data] = fluo.fluo_read_stitched_extended(fls.path.base_path, 'scanNr', scans_to_stich, 'file_type', 'raw', 'metadata_source', 'orchestra', 'dead_time_correction', true);
    spectrum_corrected = sum(sum(fluo_data_corrected, 3), 2);

    [fluo_data_raw, real_time_data, live_time_data, correction_coeff_data] = fluo.fluo_read_stitched_extended(fls.path.base_path, 'scanNr', scans_to_stich, 'file_type', 'raw', 'metadata_source', 'orchestra', 'dead_time_correction', false);
    spectrum_raw= sum(sum(fluo_data_raw, 3), 2);

    dead_time = real_time_data - live_time_data;
    dead_time = 100 * dead_time ./ real_time_data;
    spectrum_diff = spectrum_corrected ./ spectrum_raw;
    img_corrected = squeeze(sum(fluo_data_corrected(range, :, :), 1));
    img_raw = squeeze(sum(fluo_data_raw(range, :, :), 1));
    img_diff = img_corrected ./ img_raw;

    
    FigH = figure('Position', get(0, 'Screensize'));
    subplot(2, 3, 1), 
    semilogy(spectrum_raw, 'r'), 
    hold on, semilogy(spectrum_corrected, 'b'),
    legend('raw spectrum', 'corrected spectrum'),
    title('Raw and corrected spectra')

    subplot(2, 3, 2), 
    plot(spectrum_diff, 'g'), 
    title('Difference (corrected ./ raw)')

    subplot(2, 3, 3),
    imshow(dead_time), colormap bone, axis tight equal,
    caxis([min(dead_time(:)) max(dead_time(:))]), colorbar, title('Dead time map, % of real time')

    subplot(2, 3, 4),
    imshow(img_raw), colormap bone, axis tight equal,
    caxis([min(img_raw(:)) max(img_raw(:))]),
    colorbar, title('High intensity peak, not corrected')

    subplot(2, 3, 5),
    imshow(img_corrected), colormap bone, axis tight equal,
    caxis([min(img_corrected(:)) max(img_corrected(:))]),
    colorbar, title('High intensity peak, corrected')

    subplot(2, 3, 6),
    imshow(img_diff), colormap bone, axis tight equal,
    caxis([min(img_diff(:)) max(img_diff(:))]),
    colorbar, title('High intensity peak, difference (corrected ./ raw)')
    
    filename = sprintf([fls.path.output subfolder 'proj_at_theta_%03.0f.png'], sorted_theta(jj));
    saveas(FigH, filename)
    close all
end
  

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