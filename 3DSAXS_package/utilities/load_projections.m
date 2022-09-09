function [projection] = load_projections(file_name)
% The function load_projection is used to load data from sSAXS that has
% been saved with the -append command (one structure per projection)
% input: file_name: name of file saved by calling prepare_SASTT in
% saxs_caller_template in scanning_saxs
% output: the projection strcuture in the correct shape to be the input for
% SASTT analysis
% NOTE: here the q-resolved  data is ignored!

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


proj = load(file_name);
names = fieldnames(proj);
subname = fieldnames(proj.(names{1}));

projection = struct();
%reshape the input
for ii = 1:length(names)
    for iii = 1:length(subname)
        projection(ii).(subname{iii}) = proj.(names{ii}).(subname{iii});
    end
    % line added to remove nans and ifs in data
    projection(ii).data(~isfinite(projection(ii).data))=0; 
end