% function [mask3D, cut_off] = create_mask3D(tomogram, mask)
% 
% This creates a 3D binary mask based on a tomogram. Creates the map of the voxels
% to be optimized by SH
%
% Inputs:
%        - tomogram = intensity tomogram based on symmetric intensity
%        - mask = parameters to define the masking
% Output:
%        - mask3D = mask to be applied
%        - cut_off = cut off value for the intensities

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
function [mask3D, cut_off] = create_mask3D(tomogram, mask)

gauss_filter = mask.gauss_filter;
cut_off = mask.cut_off;
openning = mask.openning;
diam_cylinder = mask.diam_cylinder;

   if ~isempty(gauss_filter)
        tomogram_filt = utils.imgaussfilt3_conv(tomogram,gauss_filter);
    else
        tomogram_filt = tomogram;
    end
    if isempty(cut_off)
        [counts, x] = hist(tomogram_filt(:),200);
        cut_off = otsuthresh(counts); %values smaller than those are set to zero
    end
    fprintf('Using cut_off = %.03f \n', cut_off);
    mask3D = (tomogram > cut_off);
    
    if ~isempty(openning)
        mask3D = utils.imgaussfilt3_conv(mask3D, openning)>0.5;
    end
 
%%% Create a cylindrical mask for the sample %%%
if ~isempty(diam_cylinder)
    [ny nx nz] = size(mask3D);
    nx = [1:nx]-nx/2; % 2
    ny = [1:ny]-ny/2; % 1
    nz = [1:nz]-nz/2; % 3
    [Nx Ny Nz] = meshgrid(nx,ny,nz);
    R = sqrt(Nx.^2 + Nz.^2);
    mask3D = mask3D.*(R<diam_cylinder/2);
end  