%CONVERTER_EXIT cleanup function for convert2NEXUS
% 
% ** s              converter structure
% 
% see also: beamline.nexus.convert2NEXUS()

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


function converter_exit(s)

if isfield(s,'status')&&(s.status)
    
    fprintf('\n\n');
    
    fprintf('###############################\n');
    fprintf('###### Converter stopped ######\n');
    fprintf('###############################\n');
    
    fprintf('\n');
    
    if ~s.deletedData
        fprintf('Removing unfinished nexus file\n\t %s.\n', s.targetFile);
        try
            delete(s.targetFile);
        catch
            fprintf('Failed to remove unfinished file.\n');
        end
    end
    
    fprintf('If you want to continue, run \n');
    if s.finished
        fprintf('\t beamline.nexus.convert2NEXUS(''finished'', true, ''scanNr'', %u)\n\n', s.scanNr);
    else
        fprintf('\t beamline.nexus.convert2NEXUS(''scanNr'', %u)\n\n', s.scanNr);
    end
    
end

end