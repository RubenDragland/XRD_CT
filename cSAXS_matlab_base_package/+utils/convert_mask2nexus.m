%CONVERT_MASK2NEXUS load a .mat mask file, revert the orientation and save
%it in a nexus-like format
%
% ** sourceMask                 % full path to the .mat file
% ** targetMask                 % full path to the output .h5 file
% ** orientation                % orientation with which the .mat file was
%                               saved; should be given in [transpose fliplr flipud]
%
% EXAMPLES:
%       utils.convert_mask2nexus('./oldMask.mat', './newMask.h5', [1 1 0])
%
% see also: io.save_mask_nexus

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

function convert_mask2nexus(sourceMask, targetMask, orientation)


mask_orig = load(sourceMask);

mask = math.applyTransform(mask_orig.mask,math.convert_transform(orientation, 'transpose_fliplr_flipud', 'transpose_rot90'),'transpose_rot90',true);


io.save_mask_nexus(targetMask,mask);


end