%MANUAL_LINE_ID manually ID marked peaks
%
% ** marked_peaks         structure with peaks marked for identification
% ** energy               array with values of energies in keV identical in
%                         size to the array of channels of the detector
% ** calibration_data     data that should be used for displaying the 2D
%                         maps of the marked peaks
%
% returns:
% ++ marked_peaks         structure with additional fields with
%                         identification of the peaks
%
% EXAMPLES:
%       % mark peaks of interest
%       marked_peaks = fluo.mark_peaks_of_interest(fls);
%
%       % perform auto identification of the peaks
%       marked_peaks = fluo.manual_line_id(marked_peaks, fls.calib.energy, fluo_data);
% see also: fluo.auto_line_id()

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

function marked_peaks = manual_line_id(marked_peaks, energy, calibration_data)
got_answer = 0;
while got_answer == 0
    fsp = figure(4001);
    set(fsp, 'Name', 'Peaks identification', 'units', 'normalized', 'outerposition', [0 0.2 1 0.8])

    calibration_data_loc = calibration_data{1};
    spectrum = sum(sum(calibration_data_loc, 3), 2);

    for peak_iter = 1:numel(marked_peaks.index_center)
        clf(fsp.Number)
        subplot(2, 1, 1)
        semilogy(energy, spectrum)
        ylim([min(spectrum(:)) max(spectrum(:))*5e2])
        xlabel('Energy, [keV]')
        title(['Provide peaks identification (look into the command window).'])
        hold on
        index_center = marked_peaks.index_center(peak_iter);
        half_bandwidth = round(marked_peaks.bandwidth(peak_iter) / 2);

        mark_1_ind = index_center - half_bandwidth; % left mark
        mark_2_ind = index_center + half_bandwidth; % right mark

        H = area(energy(mark_1_ind : mark_2_ind), spectrum(mark_1_ind : mark_2_ind));
        set(H(1),'FaceColor','r');
        alpha(.3);

        value_string = sprintf('Peak #%d, E=%2.2f', peak_iter, energy(index_center));
        text_handle = text(energy(index_center), spectrum(index_center) + spectrum(index_center)/10, value_string);
        set(text_handle,'Rotation',90)

        subplot(2, 1, 2)
        image_loc = squeeze(sum(calibration_data_loc(mark_1_ind:mark_2_ind, :, :), 1));
        max_filtered = max(medfilt1(image_loc(:)));
        imshow(image_loc), axis tight equal, colormap bone, caxis([0 max_filtered])
        title(sprintf('Peak #%d, E=%2.2f', peak_iter, energy(index_center)))
        hold off

        answer = input('Please provide the name of the highlighted emission line and press Enter\nNote that this identification is not used for calibration so can be given in any form (e.g. Cu Ka or The Weird Line).\n','s');
        temp_marked_peaks_ID{peak_iter} = answer;
        temp_marked_peaks_energy{peak_iter} = energy(index_center);

    end
    marked_peaks.line_ID = {};
    for ii = 1:numel(temp_marked_peaks_ID)
        marked_peaks.line_id{ii} = temp_marked_peaks_ID{ii};
    end
    
    close(fsp.Number)
    
    
    fsp = figure(4001);
    set(fsp, 'Name', 'Confirm the identification', 'units', 'normalized', 'outerposition', [0 0 1 1])
    
    num_of_cols = ceil(sqrt(numel(marked_peaks.index_center)));
    num_of_rows = num_of_cols + 2;
    y_axis_lim = max(spectrum(:));
    
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
        
        value_string = marked_peaks.line_id{peak_iter};
        text_handle = text(energy(index_center), spectrum(index_center) + spectrum(index_center)/10, value_string);
        set(text_handle,'Rotation',90)
        hold on
    end
    
    for peak_iter = 1 : numel(marked_peaks.index_center)
        index_center = marked_peaks.index_center(peak_iter);
        half_bandwidth = round(marked_peaks.bandwidth(peak_iter) / 2);
        
        mark_1_ind = index_center - half_bandwidth; % left mark
        mark_2_ind = index_center + half_bandwidth; % right mark
        
        subplot(num_of_rows, num_of_cols, num_of_cols*2+peak_iter)
        image_loc = squeeze(sum(calibration_data_loc(mark_1_ind:mark_2_ind, :, :), 1));
        max_filtered = max(medfilt1(image_loc(:)));
        imshow(image_loc), axis tight equal, colormap bone, caxis([0 max_filtered])
        title(marked_peaks.line_id{peak_iter})
        hold on
    end

    answer = questdlg('Are you happy with your identification?', ...
        'Confirm peaks identification', ...
        'Yes, I am happy','Nope, let me ID them again','Cancel', 'Yes, I am happy');
    switch answer
        case 'Yes, I am happy'
            disp('Okay.')
            got_answer = 1;
        case 'Nope, let me ID them again'
            fprintf('\n\nSure, I am just a machine, I have plenty of time for that..\nDialing down sarcasm levels..\n')
        case 'Cancel'
            error('You have decided to quit the game without saving *SAD FACE*')
    end
    close(fsp.Number)
end

end