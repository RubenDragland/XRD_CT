function tomo_out = arb_back_projection( proj_2D, x, y, X, Y, Z, R, p)
% tomo_out = arb_back_projection( proj_2D, x, y, X, Y, Z, R, p)
%
% Computes a back projection of 2D images into 3D array after an arbitrary
% coordinate transformation defined by the rotation matrix R.
%
% Inputs
%  proj_2D        - A 2D or 3D array to allow for tensor tomography, first
%                   two indices are spatial coordinates and the third allows for instance
%                   for another quantity that varies along the fourth index, e.g. detector
%                   angular section.
%  X, Y, Z        - 3D arrays that define the array dimension that
%                   corresponds to X Y Z respectively. They must be integer
%                   valued and they must have unit interval between voxels.
%                   May be created with meshgrid. These arrays allow you to
%                   define which array dimension corresponds to X, Y, Z
%                   respectively. The rotation is around X = Y = Z = 0.
%                   Even if this point does not exist in the array.
% R               - 3x3 rotation matrix that describes the rotation of the
%                   coordinates in vector convention. [x;y;z] = R*[x;y;z].
%                   Notice this matrix may be the inverse of that used in
%                   arb_proj
% p               - Structure that defines the parameters of the projection
% p.method          'nearest' or 'bilinear', the method affects the
%                   allocation of the 2D array values.
% x, y            - Row arrays, e.g. [-10:10] that specify the input
%                   coordinates of the projection. Values should be spaced
%                   by one but a subpixel shift can be added, could be useful
%                   for including projection misalignment depending on how things
%                   are handled.
%
% Manuel Guizar - 2014.07.16

%*-------------------------------------------------------------------------------------*
%|                                                                                                           |
%|  Except where otherwise noted, this work is licensed under a            |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0        |
%|  International (CC BY-NC-SA 4.0) license.                                         |
%|                                                                                                           |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)       |
%|                                                                                                           |
%|      Author: CXS group, PSI                                                                |
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

if isfield(p,'method')
    method = p.method;
else
    method = 'bilinear';     % 'nearest' or 'bilinear' method for volume allocation
end

if (~all(size(X)==size(Y)))||(~all(size(X)==size(Z)))
    error('Size of input arrays X, Y, Z does not match')
end

% Check that spacing between voxels is 1 and that they are oriented along
% main matrix axes. This is an assumption taken for the bilinear
% interpolation
if ~((abs(X(2,2,2)-X(1,1,1))==1)&&(abs(Y(2,2,2)-Y(1,1,1))==1)&&(abs(Z(2,2,2)-Z(1,1,1))==1))
    error('Spacing between consecutive coordinates in X Y or Z is not 1')
end



% if ~((x(2)-x(1)==1)&&(y(2)-y(1)==1)) %%rounding issues
if ~(abs(y(2)-y(1))-1<1e-5)&&(abs(x(2)-x(1))-1<1e-5)
    error('x and y must be monotonically increasing in integer steps')
end

if numel(x)~=size(proj_2D,2)
    error('Size of x does not match size of proj_2D')
end
if numel(y)~=size(proj_2D,1)
    error('Size of y does not match size of proj_2D')
end


% Generate rotated coordinates
Xp = R(1,1)*X + R(1,2)*Y + R(1,3)*Z;
Yp = R(2,1)*X + R(2,2)*Y + R(2,3)*Z;
% Zp = R(3,1)*X + R(3,2)*Y + R(3,3)*Z; % Not needed


%%% Here I assign for each 4th index a tomogram and operate element by
%%% element on the 3D array. I think there should be a way to make this
%%% faster by using array operations or at least using MEX functions, but
%%% remains to be tried.

tomo_out = zeros([size(X) size(proj_2D,3)]);

for jj = 1:size(tomo_out,4)
    %%% Allocation of values on the 2D projection
    proj_2D_jj = proj_2D(:, :, jj);
    [nrows, ncols] = size(proj_2D_jj);
    
    % I did not find a way to exploit vectorial notation to make this
    % faster. Now using a loop so things may be speed up dramatically if
    % MEX functions are used.
    
    if strcmpi(method,'nearest')
        %%% Nearest neighbor %%%
        
        xi = Xp-min(x)+1;
        yi = Yp-min(y)+1;
        
        % Check for out of range values of xi and set to 1
        xout = find(xi < .5 | xi > ncols+.5);
        if ~isempty(xout)
            xi(xout) = 1;
        end
        
        % Check for out of range values of yi and set to 1
        yout = find(yi < .5 | yi > nrows+.5);
        if ~isempty(yout)
            yi(yout) = 1;
        end
        
        % Matrix element indexing
        ind = round(yi)+(round(xi)-1)*nrows;
        
        % Now interpolate
        tomo_out_aux = proj_2D_jj(ind);
        
    elseif strcmpi(method,'bilinear')
        %%% Bilinear interpolation %%%
        
        xi = Xp-min(x)+1;
        yi = Yp-min(y)+1;
        
        % Check for out of range values of xi and set to 1 (lowest possible)
        xout = find(xi < 1 | xi > ncols);
        if ~isempty(xout)
            xi(xout) = 1;
        end
        
        % Check for out of range values of yi and set to 1 (lowest possible)
        yout = find(yi < 1 | yi > nrows);
        if ~isempty(yout)
            yi(yout) = 1;
        end
        
        % Matrix element indexing
        ind = floor(yi)+(floor(xi)-1)*nrows;
        
        % Compute intepolation parameters, check for boundary value
        boundary = find(xi == ncols);
        xi = xi - floor(xi);
        if ~isempty(boundary)
            xi(boundary) = xi(boundary) + 1;
            ind(boundary) = ind(boundary) - nrows;
        end
        
        boundary = find(yi == nrows);
        yi = yi - floor(yi);
        if ~isempty(boundary)
            yi(boundary) = yi(boundary) + 1;
            ind(boundary) = ind(boundary) - 1;
        end
        
        % Now interpolate
        tomo_out_aux = (proj_2D_jj(ind).*(1-yi) + proj_2D_jj(ind+1).*yi).*(1-xi) + ...
            (proj_2D_jj(ind+nrows).*(1-yi) + proj_2D_jj(ind+nrows+1).*yi).*xi;
        
    elseif strcmpi(method, 'bicubic')
        %%% Bicubic interpolation %%%
        
        xi = Xp-min(x)+1;
        yi = Yp-min(y)+1;
        
        % Check for out of range values of xi and set to 2 (lowest possible)
        xout = find(xi < 2 | xi > ncols-1);
        if ~isempty(xout)
            xi(xout) = 2;
        end
        
        % Check for out of range values of yi and set to 2 (lowest possible)
        yout = find(yi < 2 | yi > nrows-1);
        if ~isempty(yout)
            yi(yout) = 2;
        end
        
        % Matrix element indexing
        ind = (floor(yi)-1)+(floor(xi)-2)*nrows;
        
        % Compute intepolation parameters, check for boundary value
        boundary = find(xi == ncols);
        xi = xi - floor(xi);
        if ~isempty(boundary)
            xi(boundary) = xi(boundary) + 1;
            ind(boundary) = ind(boundary) - nrows;
        end
        
        boundary = find(yi == nrows);
        yi = yi - floor(yi);
        if ~isempty(boundary)
            yi(boundary) = yi(boundary) + 1;
            ind(boundary) = ind(boundary) - 1;
        end
        
        % Now interpolate
        y0 = ((2-yi).*yi-1).*yi;
        y1 = (3*yi-5).*yi.*yi+2;
        y2 = ((4-3*yi).*yi+1).*yi;
        yi = (yi-1).*yi.*yi;
        
        tomo_out_aux = (proj_2D_jj(ind).*y0 + proj_2D_jj(ind+1).*y1 + ...
            proj_2D_jj(ind+2).*y2 + proj_2D_jj(ind+3).*yi) .* ...
            (((2-xi).*xi-1).*xi);
        ind = ind + nrows;
        tomo_out_aux = tomo_out_aux + (proj_2D_jj(ind).*y0 + proj_2D_jj(ind+1).*y1 + ...
            proj_2D_jj(ind+2).*y2 + proj_2D_jj(ind+3).*yi) .* ...
            ((3*xi-5).*xi.*xi+2);
        ind = ind + nrows;
        tomo_out_aux = tomo_out_aux + (proj_2D_jj(ind).*y0 + proj_2D_jj(ind+1).*y1 + ...
            proj_2D_jj(ind+2).*y2 + proj_2D_jj(ind+3).*yi) .* ...
            (((4-3*xi).*xi+1).*xi);
        ind = ind + nrows;
        tomo_out_aux  = tomo_out_aux + (proj_2D_jj(ind).*y0 + proj_2D_jj(ind+1).*y1 + ...
            proj_2D_jj(ind+2).*y2 + proj_2D_jj(ind+3).*yi) .* ...
            ((xi-1).*xi.*xi);
        tomo_out_aux = tomo_out_aux/4;
        
    else
        error('Method was not recognized')
    end
    
    % Set out of range values to 0
    if ~isempty(xout); tomo_out_aux(xout) = 0; end
    if ~isempty(yout); tomo_out_aux(yout) = 0; end
    
    tomo_out(:,:,:,jj) = tomo_out_aux;
end

