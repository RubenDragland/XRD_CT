% [projection] = needle_removal(projection, threshold_needle_max,threshold_needle_min , p)

%*-------------------------------------------------------------------------------------*
%|                                                                                     |
%|  Except where otherwise noted, this work is licensed under a                        |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0                          |
%|  International (CC BY-NC-SA 4.0) license.                                           |
%|                                                                                     |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)                  |
%|                                                                                     |
%|      Author: CXS group, PSI                                                         |
%*------------------------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If this code, or subfunctions or parts of it, is used for research in a
%   publication or if it is fully or partially rewritten for another
%   computing language the authors and institution should be acknowledged
%   in written form and additionally you should cite:
%     M. Liebi, M. Georgiadis, A. Menzel, P. Schneider, J. Kohlbrecher,
%     O. Bunk, and M. Guizar-Sicairos, “Nanostructure surveys of
%     macroscopic specimens by small-angle scattering tensor tomography,”
%     Nature 527, 349-352 (2015).   (doi:10.1038/nature16056)
% and
%     M. Liebi, M. Georgiadis, J. Kohlbrecher, M. Holler, J. Raabe, I.
%     Usov, A. Menzel, P. Schneider, O. Bunk and M. Guizar-Sicairos,
%     "Small-angle X-ray scattering tensor tomography: model of the
%     three-dimensional reciprocal-space map, reconstruction algorithm
%     and angular sampling requirements," Acta Cryst. A74, 12-24 (2018).
%     (doi:10.1107/S205327331701614X)
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%
% This code and subroutines are part of a continuous development, they
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [projection] = needle_removal(projection, threshold_needle_max,threshold_needle_min , p)

use_diode = p.use_diode_to_remove_needle;
se_close_rad = p.close_se_rad;
se_erode_rad = p.erode_se_rad;
whichprojections = p.projections;

if isempty(whichprojections)
    whichprojections = [1:numel(projection)];
end

if ~isempty(threshold_needle_max) || ~isempty(threshold_needle_min)
    if ~isempty(se_close_rad)
        se_close = strel('disk',se_close_rad);
    end
    if ~isempty(se_erode_rad)
        se_erode = strel('disk',se_erode_rad);
    end
    if isempty(threshold_needle_max)
        threshold_needle_max = inf; % to not exclude any values
    elseif isempty(threshold_needle_min)
        threshold_needle_min = 0; % exclude negative values, if any
    end
    for ii = [1:numel(projection)]
        if any(ii == whichprojections)
            maskneedle = ones(size(projection(ii).diode));
            projection(ii).needle_removal = [];
            if use_diode
                maskneedle(projection(ii).diode > threshold_needle_max) = 0;
                maskneedle(projection(ii).diode < threshold_needle_min) = 0;
            else
                maskneedle(mean(projection(ii).data, 3) > threshold_needle_max) = 0;
                maskneedle(mean(projection(ii).data, 3) < threshold_needle_min) = 0;
            end
            
            if ~isempty(p.max_needle_height)
                maskneedle(p.max_needle_height(ii):end,:) = 1;
            end
            if ~isempty(p.max_needle_right)
                maskneedle(:,p.max_needle_right(ii):end) = 1;
            end
            if ~isempty(p.min_needle_left)
                maskneedle(:,1:p.min_needle_left(ii)) = 1;
            end
            
            
            %   IM2 = IMERODE(IM,SE) erodes the grayscale, binary, or packed binary image
            %   IM, returning the eroded image, IM2.  SE is a structuring element
            %   object, or array of structuring element objects, returned by the
            %   STREL or OFFSETSTREL functions.
            
            if ~isempty(se_close_rad)
                maskneedle = imclose(maskneedle,se_close);
            end
            if ~isempty(se_erode_rad)
                maskneedle = imerode(maskneedle,se_erode);
            end
        else
            maskneedle = ones(size(projection(ii).diode));
        end
        % save the result
        projection(ii).needle_removal = maskneedle;
    end
    
    % plotting the results
    win_mask = [];
    for ii = 1:numel(projection)
        win_mask(:,:,ii) = projection(ii).needle_removal;
    end
    
    figure(52)
    if ~(p.use_diode_to_remove_needle)
        data_av = average_angular_sectors(projection);
        plotting.imagesc3D(data_av.*win_mask,'init_frame',numel(projection));
        colormap default;
        axis xy equal tight
    else
        diode_data = [projection.diode];
        diode_data = reshape(diode_data, size(projection(1).diode, 1), size(projection(1).diode, 2), length(projection));
        plotting.imagesc3D(diode_data.*win_mask,'init_frame',numel(projection));
        colormap default;
        axis xy equal tight
    end
    
    
else
    for ii = 1:length(projection)
        projection(ii).needle_removal = [];
    end
end