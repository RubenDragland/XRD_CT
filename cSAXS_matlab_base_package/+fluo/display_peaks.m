%DISPLAY_PEAKS function that shows 2D maps of the marked peaks, applying the 
%              marked ROIs to provided data
%
% ** fls                  fluo structure
% ** marked_peaks         structure with marked and identified peaks
% ** fluo_data            array with fluo projection (channels*y*x)
%
%
% EXAMPLES:
%       % mark peaks of interest
%       marked_peaks = fluo.mark_peaks_of_interest(fls);
%
%       % perform auto identification of the peaks
%       marked_peaks = fluo.manual_line_id(marked_peaks, fls.calib.energy, fluo_data);
%       % display marked peaks
%       scans_to_inspect = [340];
%       fluo_data = fluo.fluo_read(fls.path.base_path, 'scanNr', scans_to_inspect, 'fluo_structure', fls);
%       fluo.display_peaks(fls, marked_peaks, fluo_data)
% see also: fluo.inspect_projection()

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
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a
%   publication or if it is fully or partially rewritten for another
%   computing language the authors and institution should be acknowledged
%   in written form in the publication: “Data processing was carried out
%   using the “cSAXS software package” developed by the CXS group,
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

function display_peaks(fls, marked_peaks, fluo_data)
spectrum = fls.calib.spectrum;
energy = fls.calib.energy;
num_of_cols = ceil(sqrt(numel(marked_peaks.index_center)));
num_of_rows = num_of_cols + 2;
y_axis_lim = max(spectrum(:));

fsp = figure(4001);
set(fsp, 'Name', 'Peaks of interest.', 'units', 'normalized', 'outerposition', [0 0 1 1])
clf(fsp.Number)
subplot(num_of_rows, num_of_cols, [1:num_of_cols*2])
semilogy(energy, spectrum)
xlabel('Energy, [keV]')
ylim([0 500*y_axis_lim])
title(['Peaks of interest.'])
hold on

for peak_iter = 1 : numel(marked_peaks.index_center)
    index_center = marked_peaks.index_center(peak_iter);
    half_bandwidth = round(marked_peaks.bandwidth(peak_iter) / 2);

    mark_1_ind = index_center - half_bandwidth; % left mark
    mark_2_ind = index_center + half_bandwidth; % right mark

    H = area(energy(mark_1_ind : mark_2_ind), spectrum(mark_1_ind : mark_2_ind));
    set(H(1),'FaceColor','r');
    alpha(.3);

    value_string = sprintf('%s', marked_peaks.line_id{peak_iter});
    text_handle = text(energy(index_center), spectrum(index_center) + spectrum(index_center)/10, value_string);
    set(text_handle,'Rotation',90)
    hold on
end

fluo_data_loc = fluo_data{1};
for peak_iter = 1 : numel(marked_peaks.index_center)
    index_center = marked_peaks.index_center(peak_iter);
    half_bandwidth = round(marked_peaks.bandwidth(peak_iter) / 2);

    mark_1_ind = index_center - half_bandwidth; % left mark
    mark_2_ind = index_center + half_bandwidth; % right mark

    subplot(num_of_rows, num_of_cols, num_of_cols*2+peak_iter)
    image_loc = squeeze(sum(fluo_data_loc(mark_1_ind:mark_2_ind, :, :), 1));
    max_filtered = max(medfilt1(image_loc(:)));
    imshow(image_loc), axis tight equal, colormap bone, caxis([0 max_filtered])
    title(sprintf('%s', marked_peaks.line_id{peak_iter}))
    hold on
end

hold off
end