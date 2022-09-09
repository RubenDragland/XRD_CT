% FIND_UNUSED_GPU_CARD Go through cards until find one that is not being used
% based on memory usage
%
% gpu_id = find_unused_gpu_card
% 
% Inputs:
% ** max_card_index     Maximum index of GPU card on computer
%
% returns: 
% ++ gpu_id         Matlab ID of the unused GPU found (indices start at 1)

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

function gpu_id = find_unused_gpu_card(max_card_index)
total_memory = ones(max_card_index,1);
gpu_found = false;
gpu_id = [];
for ii = 1:max_card_index
    [users, user_memory] = utils.report_GPU_usage(ii);
    total_memory(ii) = sum(user_memory);
    if total_memory(ii) < 1400
        utils.verbose(0,sprintf('Found unused GPU card %d',ii))
        gpu_id = ii;
        gpu_found = true;
        break
    end
end

if ~gpu_found
    utils.verbose(0,'Did not find a free GPU card pausing 10 sec')
    pause(10)
end

end
