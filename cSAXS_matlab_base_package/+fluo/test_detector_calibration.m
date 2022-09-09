%TEST_DETECTOR_CALIBRATION allows you to test whether the calibration was
%                          done properly by asking you to mark peaks that
%                          you can identify, then performing auto
%                          identification from the X-ray database and then
%                          prompting you to confirm whether this
%                          identification was good
%
% ** fls             fluo structure
%
% *optional*
% ** scanNr          scan number
%
% returns:
% ++ fls             fluo structure with field 'tested' set to true
%

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

function fls = test_detector_calibration(fls, varargin)

par = inputParser;
par.addParameter('scanNr', -1, @isnumeric)

par.parse(varargin{:})
vars = par.Results;

if vars.scanNr == -1
   scanNr = fls.calib.scan_num;
else
    scanNr = vars.scanNr;
end

if ~isfield(fls.calib, 'tested')
    fls.calib.tested = false;
end

if ~fls.calib.tested
    % mark peaks for testing the calibration
    % auto ID the peaks
    % get response whether the calibration is good or not
    % if the calibration is not good then redo the calibration
    
    base_path = fls.path.base_path;
    file_type = fls.proc.file_type;
    energy = fls.calib.energy;

    fls = fluo.reduce_xray_database(fls, 'scanNr', scanNr(1));
    calibration_data = fluo.fluo_read(base_path, 'scanNr', scanNr, 'fluo_structure', fls);

    [num_of_channels, ~, ~] = size(calibration_data{1});
    channels = [1:num_of_channels];
    spectrum = zeros(num_of_channels, 1);
    for data_iter = 1 : numel(calibration_data)
        temp_data = squeeze(sum(sum(calibration_data{data_iter}, 3), 2));
        spectrum = spectrum + temp_data;
    end
    
    marked_peaks = mark_peaks(spectrum, energy, calibration_data);
    marked_peaks_auto_ID = fluo.auto_line_id(fls, marked_peaks, energy);
    display_peaks(spectrum, energy, marked_peaks_auto_ID, calibration_data);
    
    got_answer = 0;
    while got_answer == 0
    answer = input('Are you happy with the auto ID?(Y/n)\n','s');
    if strcmp(answer, 'Yes') || strcmp(answer, 'Y') || strcmp(answer, 'y') || ...
            strcmp(answer, 'yes') || strcmp(answer, '')
        marked_peaks = marked_peaks_auto_ID;
        fls.calib.tested = true;
        got_answer = 1;
    elseif strcmp(answer, 'No') || strcmp(answer, 'N') || ...
            strcmp(answer, 'n') || strcmp(answer, 'no')
        fprintf('Sorry to hear that. Then I advise you to proceed with manual ID option when marking the peaks for the reconstruction (you will find this option as you advance through the code).\n');
        fls.calib.tested = true;
        got_answer = 1;
    else
        fprintf('I do not understand this answer. Try again.\n');
        got_answer = 0;
    end
    end
end

end

function marked_peaks = mark_peaks(spectrum, energy, calibration_data)
got_answer = 0;

while got_answer == 0
    fsp = figure(4001);
    set(fsp, 'Name', 'Marking peaks', 'units', 'normalized', 'outerposition', [0 0 1 1])
    clf(fsp.Number)
    semilogy(energy, spectrum)
    xlabel('Energy, [keV]')
    
    title({'TESTING THE CALIBRATION', ...
        'Please mark peaks that YOU KNOW by clicking on the left and on the right of each peak.',...
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
    set(fsp, 'Name', 'Inspecting marked peaks', 'units', 'normalized', 'outerposition', [0 0 1 1])
    
    subplot(num_of_rows, num_of_cols, [1:num_of_cols*2])
    semilogy(energy, spectrum)
    xlabel('Energy, [keV]')
    ylim([0 20*y_axis_lim])
    title(['Peaks selected for the fine detector calibration.'])
    hold on
    
    for peak_iter = 1 : numel(marked_peaks.index_center)
        index_center = marked_peaks.index_center(peak_iter);
        half_bandwidth = round(marked_peaks.bandwidth(peak_iter) / 2);
        
        mark_1_ind = index_center - half_bandwidth; % left mark
        mark_2_ind = index_center + half_bandwidth; % right mark
        
        H = area(energy(mark_1_ind : mark_2_ind), spectrum(mark_1_ind : mark_2_ind));
        set(H(1),'FaceColor','r');
        alpha(.3);
        
        value_string = sprintf('Peak #%d', peak_iter);
        text_handle = text(energy(index_center), spectrum(index_center) + spectrum(index_center)/10, value_string);
        set(text_handle,'Rotation',90)
        hold on
    end
    
    calibration_data_loc = calibration_data{1};
    for peak_iter = 1 : numel(marked_peaks.index_center)
        index_center = marked_peaks.index_center(peak_iter);
        half_bandwidth = round(marked_peaks.bandwidth(peak_iter) / 2);
        
        mark_1_ind = index_center - half_bandwidth; % left mark
        mark_2_ind = index_center + half_bandwidth; % right mark
        
        subplot(num_of_rows, num_of_cols, num_of_cols*2+peak_iter)
        image_loc = squeeze(sum(calibration_data_loc(mark_1_ind:mark_2_ind, :, :), 1));
        max_filtered = max(medfilt1(image_loc(:)));
        imshow(image_loc), axis tight equal, colormap bone, caxis([0 max_filtered])
        title(sprintf('Peak #%d', peak_iter))
        hold on
    end
    
    answer = questdlg('Are you happy with your choice?', ...
        'Confirm peaks selection', ...
        'Ja, sehr gut','Nope, let me pick again','Cancel', 'Ja, sehr gut');
    switch answer
        case 'Ja, sehr gut'
            disp('Okay, let me see what we have gotten there..')
            got_answer = 1;
        case 'Nope, let me pick again'
            disp('Okay, try again.')
        case 'Cancel'
            error('You have decided to quit the game without saving.')
    end
    
    hold off
end
end

function display_peaks(spectrum, energy, marked_peaks, calibration_data)
num_of_cols = ceil(sqrt(numel(marked_peaks.index_center)));
num_of_rows = num_of_cols + 2;
y_axis_lim = max(spectrum(:));

fsp = figure(4001);
set(fsp, 'Name', 'Inspecting marked peaks', 'units', 'normalized', 'outerposition', [0 0 1 1])
clf(fsp.Number)

subplot(num_of_rows, num_of_cols, [1:num_of_cols*2])
semilogy(energy, spectrum)
xlabel('Energy, [keV]')
ylim([0 20*y_axis_lim])
title(['Peaks auto IDed.'])
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

calibration_data_loc = calibration_data{1};
for peak_iter = 1 : numel(marked_peaks.index_center)
    index_center = marked_peaks.index_center(peak_iter);
    half_bandwidth = round(marked_peaks.bandwidth(peak_iter) / 2);

    mark_1_ind = index_center - half_bandwidth; % left mark
    mark_2_ind = index_center + half_bandwidth; % right mark

    subplot(num_of_rows, num_of_cols, num_of_cols*2+peak_iter)
    image_loc = squeeze(sum(calibration_data_loc(mark_1_ind:mark_2_ind, :, :), 1));
    max_filtered = max(medfilt1(image_loc(:)));
    imshow(image_loc), axis tight equal, colormap bone, caxis([0 max_filtered])
    title(sprintf('%s', marked_peaks.line_id{peak_iter}))
    hold on
end

hold off
end