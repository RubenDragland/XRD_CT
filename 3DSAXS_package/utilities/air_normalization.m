% function [projection] = air_normalization(norm_air, projection, plot_fig)
% 
% This function is used when  fluctuations in flux are observed between
% projections. If the SAXS data is normalized by air counts, the mean value
% is comparable for all projection.
%
% Inputs:
%        - norm_air = number of columns to use for correction
%        - projection =  projection structure
%        - plot_fig = 1:plot the result, 0:does not plot
% Output:
%        - projection = air normalized projection structure

%*-------------------------------------------------------------------------------------*
%|                                                                                     |
%|  Except where otherwise noted, this work is licensed under a                        |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0                          |
%|  International (CC BY-NC-SA 4.0) license.                                           |
%|                                                                                     |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)                  |
%|                                                                                     |
%|      Author: CXS group, PSI                                                         |
%*------------------------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If this code, or subfunctions or parts of it, is used for research in a
%   publication or if it is fully or partially rewritten for another
%   computing language the authors and institution should be acknowledged
%   in written form and additionally you should cite:
%     M. Liebi, M. Georgiadis, A. Menzel, P. Schneider, J. Kohlbrecher,
%     O. Bunk, and M. Guizar-Sicairos, “Nanostructure surveys of
%     macroscopic specimens by small-angle scattering tensor tomography,”
%     Nature 527, 349-352 (2015).   (doi:10.1038/nature16056)
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [projection] = air_normalization(norm_air, projection, plot_fig)
Nprojections = length(projection);
for ii = 1:Nprojections
    temp = projection(ii).data;
    % select the range for air correction
    if norm_air == 1
        range = [1, size(temp, 2)];
    else
        range = [1:norm_air, size(temp, 2)-1:size(temp, 2)];
    end
    %find the median value for air
    air_edge = squeeze(median(median(temp(:, range, :))));
    if ~isfield(projection(ii).par, 'norm_air')
        projection(ii).par.norm_air =  norm_air;
    end
    % check if values change
    if projection(ii).par.norm_air ~= norm_air
        fprintf('Please reload the data to RE-normalize \n')
        return
    end
    % normalize all projections by the median air value
    for iii = 1:length(air_edge)
        temp(:,:,iii) = temp(:,:,iii)./air_edge(iii);
    end
    % replace data
    projection(ii).data = temp;
end

% Did the air normalization work?
data_ct = [];
for ii = 1:Nprojections
    data_ct = [data_ct; squeeze(mean(mean(projection(ii).data(:, :))))];
end
if plot_fig
    figure(2);
    clf
    plot(1:Nprojections, data_ct, 'ko--')
    hold on
    xlabel('projection number')
    ylabel('Mean Intensity [counts/s]')
    legend({'Normalized data'})
    grid on
    title('Did the air normalization help with SAXS data stability?')
end
