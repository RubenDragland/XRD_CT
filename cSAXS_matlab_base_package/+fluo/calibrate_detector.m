%CALIBRATE_DETECTOR calibrate detector from channels into keV
%
% ** fls             fluo structure
%
% returns:
% ++ fls             updated fluo structure with calibration values
%
% EXAMPLES:
%       fls = fluo.calibrate_detector(fls)
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

function [fls] = calibrate_detector(fls)
if ~isempty(fls.path.calibration)
    disp('NOTE: calibration file is provided, will not perform new calibration!')
    fprintf('Loading calibration data from:\n%s\n', fls.path.calibration)
    load(fls.path.calibration);
    fls.calib.scanNr = calibration.scanNr;
    fls.calib.spectrum = calibration.spectrum;
    fls.calib.energy = calibration.energy;
    fls.calib.coefs = calibration.coefs;
    fls.calib.peaks = calibration.peaks;
    
elseif isempty(fls.path.calibration) && ~isempty(fls.calib.scan_num)
    base_path = fls.path.base_path;
    scanNr = fls.calib.scan_num;
    file_type = fls.proc.file_type;

    fls = fluo.reduce_xray_database(fls, 'scanNr', scanNr(1));
    calibration_data = fluo.fluo_read(base_path, 'scanNr', scanNr, 'fluo_structure', fls);
    
    if numel(size(calibration_data{1})) == 3
        [num_of_channels, ~, ~] = size(calibration_data{1});
        channels = [1:num_of_channels];
        spectrum = zeros(num_of_channels, 1);
        for data_iter = 1 : numel(calibration_data)
            temp_data = squeeze(sum(sum(calibration_data{data_iter}, 3), 2));
            spectrum = spectrum + temp_data;
        end
    else
        error('Data should be acquired in the mapping mode (i.e. not a loopscan or alike)')
    end

    [energy, elastic_peak] = get_energy_scale_from_elastic_peak(spectrum, base_path, scanNr, file_type);
    marked_peaks = mark_peaks(spectrum, energy, calibration_data);
    [calibration_peaks, calibration_coefs] = get_calibration(fls, spectrum, energy, marked_peaks, elastic_peak, calibration_data);
    energy = calibration_coefs.P1 * channels.^2 + calibration_coefs.P2 * channels + calibration_coefs.P3;

    % exporting the calibration into substructure and also into a file
    fls.calib.scanNr = scanNr;
    fls.calib.spectrum = spectrum;
    fls.calib.energy = energy;
    fls.calib.coefs = calibration_coefs;
    fls.calib.peaks = calibration_peaks;

    file_name = fullfile([fls.path.matlab, '/base/+fluo/', 'detector_calibration_', datestr(now, 'yymmdd_HHMMSS'), '.mat']);
    calibration = fls.calib;
    fprintf('Saving detector calibration into:\n%s\n', file_name)
    save(file_name, 'calibration')
elseif isempty(fls.calib.file) && isempty(fls.calib.scan_num)
    error('Either path to calibration file or the scan numbers for calibration must be provided.')
end
end

function [energy, elastic_peak] = get_energy_scale_from_elastic_peak(spectrum, base_path, scanNr, file_type)
% get interactive input for the elastic peak for coarse calibration
fsp = figure(4001);
set(fsp, 'Name', 'Marking elastic peak', 'units', 'normalized', 'outerposition', [0 0 1 1])
clf(fsp.Number)
semilogy(spectrum)

title({'Please mark elastic peak with two points on left and right.',...
       'Pressing Backspace or Delete removes the previously selected point.',...
       'Pressing Return or Enter ends the selection.'})
[spectral_range_x, spectral_range_y] = getpts(fsp);

while (length(spectral_range_x) ~= 2)
    warndlg('Wrong number of markers. Please try again.', 'Error')
    pause(0.5);
    [spectral_range_x, spectral_range_y] = getpts(fsp);
end

close(fsp.Number)

spectral_range_x = round(spectral_range_x);
spectral_range_y = round(spectral_range_y);

[max_val, ind_max] = max(spectrum(min(spectral_range_x): max(spectral_range_x)));
index = min(spectral_range_x) + ind_max - 1;

if strcmp(file_type, 'nexus')
    beam_energy = io.nexus_read(base_path, 'scanNr', scanNr(1), 'filter', 'mokev');
elseif strcmp(file_type, 'raw')
    motor_pos = io.spec_read(base_path, 'ScanNr', scanNr(1));
    beam_energy = motor_pos.mokev;
end

energy_coeff = beam_energy / index;
energy = double(linspace(0, numel(spectrum) * energy_coeff, numel(spectrum)));

elastic_peak.energy = beam_energy;
elastic_peak.index = index;
end

function marked_peaks = mark_peaks(spectrum, energy, calibration_data)
got_answer = 0;

while got_answer == 0
    
    fsp = figure(4001);
    set(fsp, 'Name', 'Marking peaks', 'units', 'normalized', 'outerposition', [0 0 1 1])
    clf(fsp.Number)
    semilogy(energy, spectrum)
    xlabel('Energy, [keV]')
    
    title({'NOW LET''S DO THE FINE CALIBRATION', ...
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
    
    answer = questdlg('Will you be able to identify these peaks or should we try clicking again?', ...
        'Confirm peaks selection', ...
        'Yes, I can ID them','No, let me pick again','Cancel', 'Yes, I can ID them');
    switch answer
        case 'Yes, I can ID them'
            disp('Okay, let me take some paper to write down the IDs of the peaks....')
            got_answer = 1;
        case 'No, let me pick again'
            disp('Okay, try again.')
        case 'Cancel'
            error('You have decided to quit the game without saving *SAD FACE*')
    end
    
    hold off
end
end

function [calibration_peaks, calibration_coefs] = get_calibration(fls, spectrum, energy, marked_peaks, elastic_peak, calibration_data)
fsp = figure(4001);
set(fsp, 'Name', 'Peaks identification', 'units', 'normalized', 'outerposition', [0 0.2 1 0.8])

calibration_data_loc = calibration_data{1};

for peak_iter = 1:numel(marked_peaks.index_center)
    clf(fsp.Number)
    subplot(2, 1, 1)
    semilogy(energy, spectrum)
    ylim([min(spectrum(:)) max(spectrum(:))*1e2])
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
    
    value_string = sprintf('Peak #%d', peak_iter);
    text_handle = text(energy(index_center), spectrum(index_center) + spectrum(index_center)/10, value_string);
    set(text_handle,'Rotation',90)

    subplot(2, 1, 2)
    image_loc = squeeze(sum(calibration_data_loc(mark_1_ind:mark_2_ind, :, :), 1));
    max_filtered = max(medfilt1(image_loc(:)));
    imshow(image_loc), axis tight equal, colormap bone, caxis([0 max_filtered])
    title(sprintf('Peak #%d', peak_iter))
    hold off
    
    got_answer = 0;
    while got_answer == 0
        answer = input('Please provide the name of the highlighted emission line (e.g. ''Ar Ka'') and press Enter\n','s');
        answer = split(answer);
        element_name = answer(1);
        line_name = answer(2);
        
        found_lines = [];
        for database_iter = 1 : numel([fls.xray_database.reduced_lines.element])
            db_element = cellstr(fls.xray_database.reduced_lines(database_iter).element);
            db_line = cellstr(fls.xray_database.reduced_lines(database_iter).line);
            db_line = db_line{1}(1:2);
            
            if strcmp(db_element, element_name) && strcmp(db_line, line_name)
                found_lines = [found_lines database_iter];
            end
        end
        
        if ~isempty(found_lines)
            [~, max_ind] = max([fls.xray_database.reduced_lines(found_lines).intensity]);
            max_intensity_ind = found_lines(max_ind);
            line_energy = fls.xray_database.reduced_lines(max_intensity_ind).energy;
            temp_marked_peaks_energy{peak_iter} = line_energy;
            temp_marked_peaks_ID{peak_iter} = [element_name line_name];
            got_answer = 1;
        else
            fprintf('Couldn''t find the line %s %s. Please check for spelling errors.', answer(1), answer(2))
        end
    end
    
end
marked_peaks.energy = [temp_marked_peaks_energy{:}];
marked_peaks.ID = {};
for ii = 1:numel(temp_marked_peaks_ID)
    marked_peaks.ID{ii} = [temp_marked_peaks_ID{ii}{1} ' ' temp_marked_peaks_ID{ii}{2}];
end

energies_to_fit = [marked_peaks.energy elastic_peak.energy];
channels_to_fit = [marked_peaks.index_center elastic_peak.index];

fit_values = polyfit(channels_to_fit, energies_to_fit, 2);

calibration_coefs.P1 = fit_values(1);
calibration_coefs.P2 = fit_values(2);
calibration_coefs.P3 = fit_values(3);

calibration_peaks.energy = energies_to_fit;
calibration_peaks.channels = channels_to_fit;
calibration_peaks.ID = [marked_peaks.ID 'elastic peak'];

close(fsp.Number)
end