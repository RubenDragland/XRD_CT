% function [par] = find_rsum(par, fnames)
% find_rsum is used to select iteratively the srum-range to azimuthally
% integrate where the segments are strongly anisotropic
% function input
%   par = scan parameters
%   fnames = path directories
% function output
%   par = including the updated parameters to select the q-resolved

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

function [par] = find_rsum(par, fnames)

import plotting.plot_radial_integ

% use to select radio for scanning_SAXS analysis
name = fullfile(fnames.basedir_integrated_data, fnames.basename_integrated_data);
% load the data
[q, I_all] = plot_radial_integ(sprintf(name, par.test_scan), 'PlotQ', 1, 'Inverse_nm', 1, ...
    'ShowFigure', 0, 'SegAvg', 0, 'PointAvg', 1, 'XLog', 1, 'YLog', 1, 'PointRange', ...
    par.point_range, 'LegendMulSeg', 0);

% calculate mean
% load the data
[q_mean, I_mean] = plot_radial_integ(sprintf(name, par.test_scan), 'PlotQ', 1, 'Inverse_nm', 1, ...
    'ShowFigure', 0, 'SegAvg', 1, 'PointAvg', 1, 'XLog', 1, 'YLog', 1, 'PointRange', ...
    par.point_range, 'LegendMulSeg', 0);

% select q_range
fig = figure(1);
loglog(q(:,:, 1), mean(I_all, 3));
set(gca, 'XScale', 'log', 'YScale', 'log')
hold on
if (isfield(par,'PlotMulQ') && isfield(par,'QMulPow') && (par.PlotMulQ == 1))
    I_plot = q_mean.^par.QMulPow .* I_mean;
else
    I_plot = I_mean;
end
loglog(q_mean, I_plot, '.k--');
grid on
axis tight
ylabel('Intensity [a. u.]')
xlabel('q [nm^{-1}]')

% get the positions
for ii = 1:1000
    cursor = datacursormode(fig);
    set(cursor, 'DisplayStyle', 'datatip', 'SnapToDataVertex', 'off', 'Enable', 'on');
    waitforbuttonpress
    c_info = getCursorInfo(cursor);
    fprintf('q = %0.2d nm^{-1},  d = %0.2f nm,   r_sum = %d \n', c_info.Position(1), 6.28/c_info.Position(1), c_info.DataIndex)
end
