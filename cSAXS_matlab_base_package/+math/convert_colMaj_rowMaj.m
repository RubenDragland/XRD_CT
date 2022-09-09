%CONVERT_COLMAJ_ROWMAJ Convert integration mask indices from column-major
%to row-major
% 
% ** integ_masks            integration mask structure
% ** asize                  array dimensions
%
% returns:
% ++ integ_masks            updated integ_masks structure with indices in
%                           row-major order
%
% NOTE: The function can also be used to convert from row-major to
% column-major by giving a row-major integration mask structure as input
% with a transposed asize.
%
% see also: beamline.convert_integ_mask

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


function integ_masks = convert_colMaj_rowMaj(integ_masks, asize)

radii = size(integ_masks.indices,1);
segs = size(integ_masks.indices,2);

dims = asize;

for seg=1:segs
    for rad=1:radii
        tmpIndices = integ_masks.indices{rad,seg};
        for ii =1:length(tmpIndices)
            if tmpIndices(ii) > 0
                tmpIndices(ii) = mod(tmpIndices(ii)-1, dims(1))*dims(2)+ceil(tmpIndices(ii)/dims(1));
            end
        end
        integ_masks.indices{rad,seg} = tmpIndices;
    end
end

end