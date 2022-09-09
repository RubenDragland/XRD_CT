%% UPDATE_MASK
% This small script guides you to update an alread existing mask for
% ptychography. The main tool for creating a new mask is
% beamline.create_mask, a GUI that lets you select bad/hot pixels. 
% UPDATE_MASK loads the data, specified by file_path, plots it and starts
% the GUI. Although you can create a 3D mask, i.e. a mask which varies from
% frame to frame, a 2D mask is sufficient for most datasets.

% You can load an already existing mask within the GUI.

close all

nxs = io.nexus_read('~/Data10', 'scanNr', 98, 'filter', 'eiger_4','orientation',[0 0]);
nxs_dark = io.nexus_read('~/Data10', 'scanNr', 78, 'filter', 'eiger_4','orientation',[0 0]);


%% plot the data
figure(1),
plotting.imagesc3D(abs(log10(double(nxs.data)+1)));
colorbar
axis xy equal tight
colorbar 
title('Detector raw data')
colormap(plotting.colormaps.fire_map)

%% Optional - create a mask based on dark frames
% Here you create a mask in a temporary file, then you can load it
% additionally to your default mask in the next step
create_dark_frame_mask = true;

if create_dark_frame_mask
    mask_filename = 'temp_mask.mat';
    mask = any(nxs_dark.data < 1,3); %#ok<NASGU> % If we put < 2 as the threshold it captures many cosmic rays
    save(mask_filename,'mask')
    fprintf('Saving temporary mask: %s\n',mask_filename);
end
   

%% iterative step for a mask update (add dead pixels to the current mask)
mask = beamline.create_mask;


%% check it again
figure (2),
imagesc(mask); axis equal tight xy
title('Final mask')

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)    |
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

