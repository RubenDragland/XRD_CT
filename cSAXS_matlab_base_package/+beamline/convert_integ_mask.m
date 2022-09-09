%CONVERT_INTEG_MASK convert integration mask from cell to maskPointer/maskIndex pair
% 
% ** integ_masks        structure containing the integration mask; either
%                       as cell or as maskPointer/maskIndex pair
%
% returns:
% ++ integ_masks        updated integ_masks structure
%
% NOTE: The maskPointer/value pair is designed to improve compatibility with
% other programming languages and minimize disk space. "maskIndex" contains 
% all indices of the integration mask. "maskPointer" is a 3D
% array with dimensions (2, segments, radii). For reach radius and segment
% the entry defines the position of the first index within maskIndex and 
% the number of values to read. 
% Reading the values for radius "radius_ii" and segment "segm_jj" can
% therefore be done using:
%   res = maskIndex(maskPointer(radius_ii,segm_jj,1)+1:maskPointer(radius_ii,segm_jj,1)+maskPointer(radius_ii,segm_jj,2));
%
% see also: beamline.prep_integ_masks_nexus.m

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

function integ_masks = convert_integ_mask(integ_masks)

if size(integ_masks.indices,1)>1
    toCell = false;
else
    toCell = true;
end

useWeights = false;
if toCell
    if isfield(integ_masks.indices, 'maskWeight')
        useWeights = true;
    end
else
    if isfield(integ_masks, 'weights')
        useWeights = true;
    end
end

if ~toCell
    % convert from cell to maskPointer/maskIndex
    radii = size(integ_masks.indices,1);
    segments = size(integ_masks.indices,2);
    maskPointer = zeros(radii,segments,2);
    weights = [];
    flatMask = [];
    indx = 0;
    for ii=1:radii
        for jj=1:segments
            maskPointer(ii,jj,:) = [indx length(integ_masks.indices{ii,jj})];
            indx = indx + length(integ_masks.indices{ii,jj});
            if ~isempty(integ_masks.indices{ii,jj})
                flatMask = [flatMask integ_masks.indices{ii,jj}(:)'];
                if useWeights
                    weights = [weights integ_masks.weights{ii,jj}(:)'];
                end
            end
        end
    end
    integ_masks.indices = [];
    integ_masks.indices.maskPointer = permute(maskPointer, [3 2 1]);
    integ_masks.indices.maskIndex = flatMask;
    if useWeights
        integ_masks.indices.maskWeight = weights;
    end

else
    % convert back to a Matlab cell
    maskPointer = permute(integ_masks.indices.maskPointer, [3 2 1]);
    radii = size(maskPointer,1);
    segments = size(maskPointer,2);
    indices = cell(radii, segments);
    weights = cell(radii, segments);
    for ii=1:radii
        for jj=1:segments
            if maskPointer(ii,jj,2)>0
                indices{ii,jj} = integ_masks.indices.maskIndex(maskPointer(ii,jj,1)+1:maskPointer(ii,jj,1)+maskPointer(ii,jj,2));
                if useWeights
                    weights{ii,jj} = integ_masks.indices.maskWeight(maskPointer(ii,jj,1)+1:maskPointer(ii,jj,1)+maskPointer(ii,jj,2));
                end
            end
        end
    end
    integ_masks.indices = [];
    integ_masks.indices = indices;
    if useWeights
        integ_masks.weights = weights;
    end
end

end
