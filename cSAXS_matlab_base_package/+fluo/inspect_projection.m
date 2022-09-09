%INSPECT_PROJECTION function that allows you to mark peaks and inspect
%                   their 2D maps in provided dataset
%
% ** fls                  fluo structure
% ** fluo_data            array with fluo projection (channels*y*x)
%
% EXAMPLES:
%       fluo_data = fluo.fluo_read(fls.path.base_path, 'scanNr', 340, 'fluo_structure', fls);
%       fluo.inspect_projection(fluo_data, fls)
% see also: fluo.display_peaks()

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

function inspect_projection(fluo_data, fls)
fluo_data = fluo_data{1};
spectrum = sum(sum(fluo_data, 3), 2);
if isfield(fls.calib,'energy')
    energy = double(fls.calib.energy);
    channels = false;
else
    energy = [1:numel(spectrum)];
    channels = true;
end

got_answer = 0;

while got_answer == 0
    fsp = figure(4001);
    set(fsp, 'Name', 'Marking peaks', 'units', 'normalized', 'outerposition', [0 0 1 1])
    clf(fsp.Number)
    semilogy(energy, spectrum)
    if ~channels
        xlabel('Energy, [keV]')
    else
        xlabel('Channel number, [ ]')
    end
    
    title({'INSPECTING DATA', ...
        'Please mark peaks by clicking on the left and on the right of each peak.',...
        'Pressing Backspace or Delete removes the previously selected point.',...
        'Pressing Return or Enter ends the selection.'})
    [spectral_range_x, spectral_range_y] = getpts(fsp);
    
    while mod(numel(spectral_range_x),2)
        warndlg('Number of markers needs to be even. Please try again.', 'Error')
        pause(0.5);
        [spectral_range_x, spectral_range_y] = getpts(fsp);
    end
    close(fsp.Number)
    
    marked_peaks = [];
    for ch_iter = 1:length(spectral_range_x)/2
        mark_1 = spectral_range_x((ch_iter - 1)*2 + 1);
        mark_2 = spectral_range_x((ch_iter - 1)*2 + 2);
        marks = [mark_1 mark_2];
        
        mark_1 = min(marks); % left mark
        mark_2 = max(marks); % right mark
        
        ind_range = find(energy > mark_1 & energy < mark_2);
        
        [~, ind_max] = max(spectrum(ind_range));
        index_center = min(ind_range) + ind_max - 1;
        bandwidth = max(ind_range) - min(ind_range);
        marked_peaks.index_center{ch_iter} = index_center;
        marked_peaks.bandwidth{ch_iter} = bandwidth;
    end
    
    [~, indices] = sort([marked_peaks.index_center{:}]);
    marked_peaks.index_center = [marked_peaks.index_center{indices}];
    marked_peaks.bandwidth = [marked_peaks.bandwidth{indices}];
    
    num_of_cols = ceil(sqrt(numel(marked_peaks.index_center)));
    num_of_rows = num_of_cols + 2;
    y_axis_lim = max(spectrum(:));
    
    fsp = figure(4001);
    set(fsp, 'Name', 'Inspecting the peaks', 'units', 'normalized', 'outerposition', [0 0 1 1])
    clf(fsp.Number)
    subplot(num_of_rows, num_of_cols, [1:num_of_cols*2])
    semilogy(energy, spectrum)
    xlabel('Energy, [keV]')
    ylim([0 500*y_axis_lim])
    title(['Peaks selected for the inspection.'])
    hold on
    
    for peak_iter = 1 : numel(marked_peaks.index_center)
        index_center = marked_peaks.index_center(peak_iter);
        half_bandwidth = round(marked_peaks.bandwidth(peak_iter) / 2);
        
        mark_1_ind = index_center - half_bandwidth; % left mark
        mark_2_ind = index_center + half_bandwidth; % right mark
        
        H = area(energy(mark_1_ind : mark_2_ind), spectrum(mark_1_ind : mark_2_ind));
        set(H(1),'FaceColor','r');
        alpha(.3);
        
        value_string = sprintf('Peak #%d, E=%2.2fkeV', peak_iter, energy(index_center));
        text_handle = text(energy(index_center), spectrum(index_center) + spectrum(index_center)/10, value_string);
        set(text_handle,'Rotation', 90)
        hold on
    end
    
    calibration_data_loc = fluo_data;
    for peak_iter = 1 : numel(marked_peaks.index_center)
        index_center = marked_peaks.index_center(peak_iter);
        half_bandwidth = round(marked_peaks.bandwidth(peak_iter) / 2);
        
        mark_1_ind = index_center - half_bandwidth; % left mark
        mark_2_ind = index_center + half_bandwidth; % right mark
        
        subplot(num_of_rows, num_of_cols, num_of_cols*2+peak_iter)
        image_loc = squeeze(sum(calibration_data_loc(mark_1_ind:mark_2_ind, :, :), 1));
        max_filtered = max(medfilt1(image_loc(:)));
        imshow(image_loc), axis tight equal, colormap bone, caxis([0 max_filtered])
        title(sprintf('Peak #%d, E=%2.2fkeV', peak_iter, energy(index_center)))
        hold on
    end
    
    answer = questdlg('Are you happy with the selection?', ...
        'Confirm peaks selection', ...
        'Yes','No, let me pick again','Cancel', 'Yes');
    switch answer
        case 'Yes'
            got_answer = 1;
        case 'No, let me pick again'
            disp('Okay, try again.')
        case 'Cancel'
            error('You have decided to quit the game without saving *SAD FACE*')
    end
    
    hold off
end

end