%FIND_NEXT_NEXUS_FILE find next nexus file
% 
% ** scanNr             current scan number
%
% *optional*
% ** dataDir            name of the nexus directory; default: 'data'
% ** basePath           path tot the nexus directory; default: '~/Data10';
%
% returns:
% ++ fpath              path to the nexus file
% ++ nextScanNr         next scan number
%
% EXAMPLES:
%       fpath = utils.find_next_nexus_file(25);  
%       [fpath, nextScanNr] = utils.find_next_nexus_file(25, 'basePath','/das/work/p16/p16812/online');
%
%
% see also: io.nexus_read()

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


function [fpath, nextScanNr] = find_next_nexus_file(scanNr, varargin)

par = inputParser;
par.addParameter('dataDir', 'data', @ischar)
par.addParameter('basePath', '~/Data10', @ischar)

par.parse(varargin{:})
vars = par.Results;

nextScanNr = [];
resParentID = [];

% check parent directory (e.g. S01000-01999)
resParent = dir(fullfile(vars.basePath, vars.dataDir));
for ii=1:length(resParent)
   if strcmpi(utils.compile_x12sa_dirname(scanNr, true), [resParent(ii).name '/']) 
       resParentID = ii;
   end
end

if isempty(resParentID)
    % given scan number does not fit in any of these parent directories
    error('Could not find scan S%05d.', scanNr);
end


resSub = dir(fullfile(vars.basePath, vars.dataDir, utils.compile_x12sa_dirname(scanNr,true)));
scanID = 1;
while scanID <= length(resSub)
    if length(resSub(ii).name)==6
        currScanNr = str2double(resSub(scanID).name(2:end));
        if currScanNr > scanNr
            nextScanNr = currScanNr;
            break;
        end
    end
    scanID = scanID + 1;
end


% check next parent directory
if isempty(nextScanNr) && (resParentID < length(resParent))
    res = dir(fullfile(vars.basePath, vars.dataDir, resParent(resParentID+1).name));
    if ~isempty(res)
        scanID = 1;
        while scanID <= length(res)
            if length(res(scanID).name)==6
                nextScanNr = str2double(res(scanID).name(2:end));
                break;
            end
            scanID = scanID + 1;
        end
    end
end


if ~isempty(nextScanNr)
    fname = dir(fullfile(vars.basePath, vars.dataDir, utils.compile_x12sa_dirname(nextScanNr), '*.h5'));
    fpath = fullfile(fname(1).folder, fname(1).name);
else
    fpath = [];
end


    

end