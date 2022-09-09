% function [par] = find_qrange(par, fnames, plot_qres)
% find_qrange is used to select a range of q-vectors to perform q-resolved
% SASTT on
% function input
%   par = scan parameters
%   fnames = path directories
%   plot_qres = 0: if not plotting is requested, = 1: if plotting is
%   requested
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

function [par] = find_qrange(par, fnames, plot_qres)

import plotting.plot_radial_integ

% check number of radius in q-resolved
if ((par.pixel_range(end)-par.pixel_range(1))/par.int_integ) < 3
    warning('There are less than 3 pixel radii in the whole q-range selected (par.pixel_range). Consider decreasing int_integ value')
end

% check number of radius in q-resolved
if ((par.pixel_range(end)-par.pixel_range(1))/par.int_integ) < 1
    error('There is less than 1 pixel radii in the whole q-range selected (par.pixel_range). Decrease int_integ value')
end

% load the radial integration
name = fullfile(fnames.basedir_integrated_data, fnames.basename_integrated_data);
% load the data
[~, I_all] = plot_radial_integ(sprintf(name, par.test_scan), 'PlotQ', 1, 'Inverse_nm', 1, ...
    'ShowFigure', 0, 'SegAvg', 1, 'PointAvg', 1, 'XLog', 1, 'YLog', 1, 'PointRange', ...
    par.point_range, 'LegendMulSeg', 0);

% get values for q-resolved
if (par.log_scale)
    qq = round(logspace(log10(par.pixel_range(1)),log10(par.pixel_range(end)), par.int_integ));
else
    qq = round(linspace((par.pixel_range(1)),(par.pixel_range(end)), par.int_integ));
end
range_q = [par.pixel_range(1):par.pixel_range(end)];

% to replace with new data
if isfield(par, 'qresolved_q')
    warning('Field qresolved_q exists! This will replace existing values')
    par = rmfield(par, 'qresolved_q');
end

for ii = 1:numel(qq)-1
    % find its equivalent in q and provide all points in the q range
    par.qresolved_q{ii} = [range_q(range_q == qq(ii)):range_q(range_q == qq(ii+1))];
end

%%% The following is to check and remove potentially repeated q-range bins
aux = par.qresolved_q;
par.qresolved_q = {};
counter = 2;
par.qresolved_q{1} = aux{1};
for ii = 2:numel(aux)
    if ( numel(aux{ii})==numel(aux{ii-1}) ) && (all(aux{ii} == aux{ii-1}))
        % they are the same, skip them
    else
        par.qresolved_q{counter} = aux{ii};
        counter = counter+ 1;
    end
end

% plot
if plot_qres == 1
    fig = figure(1);
    clf(fig)
    loglog(range_q, I_all(range_q), 'o--r');
    set(gca, 'XScale', 'log', 'YScale', 'log')
    hold on
    loglog(qq, mean(I_all(range_q)), '.k')
    axis tight
    grid on
    ylabel('Intensity [a. u.]')
    xlabel('Radius')
    title(sprintf('%d Intervals', par.int_integ))
end