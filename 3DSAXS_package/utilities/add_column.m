% function [projection, data1] = add_column(projection, data1)
% 
% This adds a column to SAXS, transmissiona and needle mask data
% if the number of columns  are uneven. This is still needed for ASTRA alignment
%
% Inputs:
%        - projection = projection structure
%        - data1 =  input data for alignment
% Output:
%        - projection = corrected projection structure
%        - data1 =  corrected input data for alignment

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
function [projection, data1] = add_column(projection, data1)

% get the x_dimension to the even size: symmetry condition for designFilter
if mod(size(projection(1).data, 2),2) ~= 0
    % if number is odd, repeat last column in X
    if exist('data1')
        data1 = [data1 data1(:,end,:)];
    end
    Nprojections = length(projection);
    for ii = 1:Nprojections
        projection(ii).data = [projection(ii).data projection(ii).data(:,end,:)];
        if isfield(projection(1),'diode')
            projection(ii).diode = [projection(ii).diode projection(ii).diode(:,end,:)];
        end
        if isfield(projection(1),'needle_removal')&&(~isempty(projection(1).needle_removal))
            if mod(size(projection(ii).needle_removal, 2),2) ~= 0
                projection(ii).needle_removal = [projection(ii).needle_removal, projection(ii).needle_removal(:,end)];
            end
        end
    end
end
end

