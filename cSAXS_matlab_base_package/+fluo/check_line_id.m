%CHECK_LINE_ID function that checks whether the auto identification was
%succesfull. The functions shows the maps of particular lines and asks user
%whether the ID is correct.
%
% ** fls                  fluo structure
% ** marked_peaks         structure with marked and identified peaks
% ** energy               array with values of energies in keV identical in
%                         size to the array of channels of the detector
% ** spectrum             array of intensity values
%
% returns:
% ++ marked_peaks         structure where identification field is changed
%                         to 'No ID' for the lines that user labeled as
%                         misidentified
%
% EXAMPLES:
%       % mark peaks of interest
%       marked_peaks = fluo.mark_peaks_of_interest(fls);
%
%       % perform auto identification of the peaks
%       marked_peaks = fluo.manual_line_id(marked_peaks, fls.calib.energy, fluo_data);
%       % check identification
%       marked_peaks = fluo.check_line_id(fls, marked_peaks, fls.calib.spectrum, fls.calib.energy)
% see also: fluo.auto_line_id(), fluo.manual_line_id()

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

function [marked_peaks] = check_line_id(fls, marked_peaks, spectrum, energy)
fsp = figure(4001);

for energy_iter = 1:numel(marked_peaks)
    clf
    semilogy(energy, spectrum)
    xlabel('Energy, [keV]')
    title(['Check auto peak indentification.'])
    hold on
    index_center = marked_peaks{energy_iter}.index_center;
    half_bandwidth = round(marked_peaks{energy_iter}.bandwidth / 2);
    
    mark_1_ind = index_center - half_bandwidth; % left mark
    mark_2_ind = index_center + half_bandwidth; % right mark
    
    H = area(energy(mark_1_ind : mark_2_ind), spectrum(mark_1_ind : mark_2_ind));
    set(H(1),'FaceColor','r');
    alpha(.3);
    
    value_string = sprintf('%s', marked_peaks{energy_iter}.line_id);
    text(energy(index_center), spectrum(index_center) + spectrum(index_center)/10, value_string);
    hold off
    
    got_answer = 0;
    while got_answer == 0
    answer = input('is the line ID correct?(Y/n)\n','s');
    if strcmp(answer, 'Yes') || strcmp(answer, 'Y') || strcmp(answer, 'y') || ...
            strcmp(answer, 'yes') || strcmp(answer, '')
        got_answer = 1;
    elseif strcmp(answer, 'No') || strcmp(answer, 'N') || ...
            strcmp(answer, 'n') || strcmp(answer, 'no')
        marked_peaks{energy_iter}.line_id = 'No ID';
        marked_peaks{energy_iter}.xray_database_index = [];
        got_answer = 1;
    else
        fprintf('I do not understand this answer. Try again.\n');
        got_answer = 0;
    end
    end
end
close(fsp)
end