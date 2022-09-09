% function [ data ] = average_angular_sectors( projection )
% Receives a set of projections prepared for SASTT and computes the average
% of the sectors, i.e. the total scattering, using the normalized sum
% factor to account for the number of pixels in each sector.
%   Important note, it will assume that the integration masks are the same
%   and only compute the normalization factor using the mask of the first
%   projection.
%
% Input
%   projection  Stack of SASTT projections
% Output
%   data        3D array with the angular sector average of all
%               projections.

%*------------------------------------------------------------------------*
%|                                                                        |
%|  Except where otherwise noted, this work is licensed under a           |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0             |
%|  International (CC BY-NC-SA 4.0) license.                              |
%|                                                                        |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)     |
%|                                                                        |
%|      Author: CXS group, PSI                                            |
%*------------------------------------------------------------------------*
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ data ] = average_angular_sectors( projection )

size_data = maxsize(projection);
if mod(size_data(2),2)>0
    size_data(2) = size_data(2)+1;
end
data = zeros([size_data(1:2) numel(projection)]);
norm_sum = projection(1).integ.norm_sum_avg;
norm_sum = repmat(sum(norm_sum, 2), 1, size_data(1), size_data(2))./sum(norm_sum(:));
norm_sum = permute(norm_sum, [2 3 1]);
for ii = 1:numel(projection)
    dsize = size(projection(ii).data);
    norm_sum = projection(ii).integ.norm_sum_avg;
    norm_sum = repmat(sum(norm_sum, 2), 1, dsize(1), dsize(2))./sum(norm_sum(:));
    norm_sum = permute(norm_sum, [2 3 1]);
    data(1:dsize(1),1:dsize(2),ii) = sum(projection(ii).data.*norm_sum, 3);
end

end

