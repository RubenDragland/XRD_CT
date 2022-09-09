% [projection] = window_mask(projection, threshold_valid, mask_horiz, mask_vert)

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

function [projection] = window_mask(projection, threshold_valid, mask_horiz, mask_vert)

pp.volume_upsampling = 1;  % Upsamples the volume by a factor of 2 to reduce artifacts
pp.method = 'bilinear';    % 'nearest' or 'bilinear' method for volume allocation
pp.filter_2D = 3;          % Strength of filter applied to image

% create a mask for the projection, this should be the same size as the
% tomogram in the end
a_mask = ones(size(projection(1).diode, 1),size(projection(1).diode, 2),size(projection(1).diode, 2));
NN = size(a_mask);
% create the mesh for arb_projection
xx = [1:NN(2)]-ceil(NN(2)/2);
yy = [1:NN(1)]-ceil(NN(1)/2); % invert X and Y
zz = [1:NN(3)]-ceil(NN(3)/2);
[XX, YY, ZZ] = meshgrid(xx,yy,zz); % (x,y,z) (2,1,3)

for ii = 1:length(projection)
    fprintf('Calculating window mask: projection %d/%d \n', ii, length(projection))
    projection(ii).window_mask= [];
    projection(ii).window_mask = projection(ii).needle_removal;
    %load the rotation matrix
    
    R = projection(ii).Rot_exp;
    
    %%% Generate a volume and coordinates for each projection to allow
    %%% different FOV for the scanning, probably this should be a function
    % include the delta from alignment
    xxout = [1:size(projection(ii).diode,2)] - ...
        ceil(size(projection(ii).diode,2)/2);
    yyout = [1:size(projection(ii).diode,1)] - ...
        ceil(size(projection(ii).diode,1)/2);
    
    % recalculate the projection of the mask and replace values
    [proj_mask] = arb_projection(a_mask, XX, YY, ZZ, R, pp, xxout, yyout);
    
    % new mask including the values of the registration
    mask = ones(length(yyout),length(xxout));
    % make a window (less than 20)
    mask(proj_mask < threshold_valid) = 0;
%     figure(100); imagesc(proj_mask); axis xy equal tight; drawnow;
    % y, x
    if mask_horiz
        mask(:,end-mask_horiz+1:end)=0;
        mask(:,1:mask_horiz)=0;
    end
    if mask_vert
        mask(end-mask_vert+1:end,:)=0;
        mask(1:mask_vert,:)=0;
    end
    
    % for negative data values, put mask value to 0
    if isfield(projection(ii),'window_mask')&&~isempty(projection(ii).window_mask)
        projection(ii).window_mask = projection(ii).window_mask .* mask;
    else
        projection(ii).window_mask = mask;
    end
    projection(ii).window_mask(projection(ii).data < 0) = 0;
    
end