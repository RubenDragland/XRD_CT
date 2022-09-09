% APPLY_3D_APODIZATION Smoothly apodize tomogram to avoid sharp edges and air affecting
% the FRC analysis 
%
% [tomogram,circulo] = apply_3D_apodization(tomogram, rad_apod, axial_apod, radial_smooth)
%
% Inputs:
%     **tomogram - volume to be apodized 
%     **rad_apod - number of pixels to be zeroed from edge of the tomogram 
%     **axial_apod - roughly number of pixels to be zeroed from top / bottom 
%     **radial_smooth - smoothness of the apodization in pixels, default = Npix/10
% Outputs: 
%     ++tomogram - apodized volume 
%     ++circulo -apodization mask 

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


function [tomogram,circulo] = apply_3D_apodization(tomogram, rad_apod, axial_apod, radial_smooth)
    import utils.*
    [Npix,~,Nlayers] = size(tomogram);

    if nargin < 4 
        radial_smooth = Npix/10;
    end
    
    if nargin < 3
        axial_apod = [];
    end
    if ~isempty(rad_apod) 
        xt = -Npix/2:Npix/2-1;
        [X,Y] = meshgrid(xt,xt);
        radial_smooth = max(radial_smooth,1); % prevent division by zero 
        circulo= single(1-radtap(X,Y,radial_smooth,round(Npix/2-rad_apod-radial_smooth)));  
        tomogram = bsxfun(@times, tomogram, circulo);
    end
    if ~isempty(axial_apod) && Nlayers > 1
        filters = fract_hanning_pad(Nlayers,Nlayers,max(0,round(Nlayers-2*axial_apod)));
        filters = ifftshift(filters(:,1));
        tomogram = bsxfun(@times,tomogram,reshape(filters,1,1,[]));
    end
end