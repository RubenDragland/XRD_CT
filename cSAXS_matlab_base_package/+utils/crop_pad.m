% CROP_PAD adjusts the size by zero padding or cropping
% Inputs: 
%   **img                input image
%   **outsize            size of final image
% *optional:*
%   **fill               value to fill padded regions 
% returns: 
%   ++imout              cropped image 

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

function [ imout ] = crop_pad( img, outsize, fill)

if nargin < 1
    fprintf('CROP_PAD: adjusts the size by zero padding or cropping\n');
    fprintf('crop_pad(img, outsize)\n');
    return
end

Nin = size(img);

if isempty(outsize) || all(outsize(1:2) == Nin(1:2))
    imout = img;   % if outsize == [], return the same image without changes 
    return
end

Nout = outsize(1:2);

if nargin < 3
    fill = 0;
end



center = floor(Nin(1:2)/2)+1;

imout = zeros([Nout,Nin(3:end)],'like',img);

if fill ~= 0 
    imout = imout + fill; 
end

centerout = floor(Nout/2)+1;

cenout_cen = centerout - center;
imout(max(cenout_cen(1)+1,1):min(cenout_cen(1)+Nin(1),Nout(1)),max(cenout_cen(2)+1,1):min(cenout_cen(2)+Nin(2),Nout(2)),:,:) ...
    = img(max(-cenout_cen(1)+1,1):min(-cenout_cen(1)+Nout(1),Nin(1)),max(-cenout_cen(2)+1,1):min(-cenout_cen(2)+Nout(2),Nin(2)),:,:);

if ~isreal(img)
    imout = complex(imout);
end

end

