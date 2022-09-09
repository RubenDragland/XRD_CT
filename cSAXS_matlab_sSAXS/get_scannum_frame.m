% GET_SCANNUM_FRAME  Get scan numbers and spec metadata of first scan from frame ID
%
% [scann, p] = get_scannum_frame(base_path,frameid)
% 
% Inputs: 
%  ** base_path    e.g.  '~/Data10'
%  ** frameid      Frame ID number
% 
% *returns*
%  ++ scann     Scan numbers corresponding to frameid
%  ++ p         Spec metadata of the first scan in scann
%
% EXAMPLE:
%
% [scann, p] = get_scannum_frame('~/Data10', 3)

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
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
%   using the “cSAXS matlab package” developed by the CXS group,
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

function [scann, p] = get_scannum_frame(base_path,frameid)
    strToFind = sprintf('#C meta frame_id int 1 %d',frameid);
    scann = utils.spec_find_scans(base_path,strToFind);
    if (~isempty(scann))
        p = io.spec_read(base_path,'ScanNr',scann(1));
    else
        p = [];
    end
end
