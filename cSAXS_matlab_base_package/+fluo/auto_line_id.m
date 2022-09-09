%AUTO_LINE_ID automatically ID marked peaks
%
% ** fls                  fluo structure
% ** marked_peaks         structure with peaks marked for identification
% ** energy               array with values of energies in keV identical in
%                         size to the array of channels of the detector
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
%       marked_peaks = fluo.auto_line_id(fls, marked_peaks, fls.calib.energy);
% see also: fluo.manual_line_id()

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

function marked_peaks = auto_line_id(fls, marked_peaks, energy)

for peaks_iter = 1 : numel(marked_peaks.index_center)
    peak_energy = energy(marked_peaks.index_center(peaks_iter));
    [line_id, xray_database_index] = get_line_id(fls, peak_energy);
    
    marked_peaks.line_id{peaks_iter} = line_id;
    marked_peaks.xray_database_index(peaks_iter) = xray_database_index;
end

end


function [line_id, xray_database_index] = get_line_id(fls, peak_energy)
energies = [fls.xray_database.reduced_lines.energy];
index = find(abs(energies - peak_energy) == min(abs(energies - peak_energy)), 1);
line_id = sprintf('%s %s %2.2fkeV', fls.xray_database.reduced_lines(index).element,...
                                    fls.xray_database.reduced_lines(index).line,...
                                    fls.xray_database.reduced_lines(index).energy);
xray_database_index = index;
end