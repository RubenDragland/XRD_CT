%UNPACK_XRAY_DATABASE    function to convert X-ray database from .mat file
%                        where it is kept as table to structure. Structure
%                        needs more memory to keep in .mat file so this
%                        function is more optimal for repository size.
%
% ** fls             fluo structure with path to X-ray database
%
% returns:
% ++ fls             fluo structure with unpacked database
%

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

function [fls] = unpack_xray_database(fls)
[edges, lines] = extract_structure(fls);
fls.xray_database.edges = edges;
fls.xray_database.lines = lines;
end

function [edges, lines] = extract_structure(fls)
load(fls.path.xray_database);
edges = table2struct(edges);
lines = table2struct(lines);
end