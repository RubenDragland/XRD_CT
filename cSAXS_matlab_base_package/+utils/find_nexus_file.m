%FIND_NEXUS_FILE find absolute path to nexus file
% 
% ** fpath              parent directory of a nexus file
% ** scanNr             scan number
%
% returns:
% ++ fpath              path to the nexus file
%
% EXAMPLE:
%       fpath = utils.find_nexus_file('~/Data10', 25);  
%
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

function fpath = find_nexus_file(fpath, varargin)
    
    if nargin>1
        scanNr = varargin{1};
    else
        scanNr = -1;
    end
    
    if exist(fpath, 'file')==2
        return 
    else
        [~, ~, fparts] = fileparts(fpath);
        if ~isempty(fparts)
            error('File does not exist.')
        end
        if scanNr < 0
            error('Please specify a valid scan number.')
        end
    end
    
    x12saDir = utils.compile_x12sa_dirname(scanNr);
    x12saDir = strsplit(x12saDir, '/');
    
    while exist(fpath, 'file')~=2
        % check current directory
        nexusFiles = dir(fullfile(fpath, sprintf('*%u*.h5', scanNr)));
        if isempty(nexusFiles)
            mod = false; % keep track of modifications
            res = dir(fpath);
            for ii=1:length(res)
                switch res(ii).name
                    case 'Data10'
                        fpath = fullfile(fpath, 'Data10');
                        mod = true;
                    case 'data'
                        fpath = fullfile(fpath, 'data');
                        mod = true;
                    case x12saDir{1}
                        fpath = fullfile(fpath, x12saDir{1});
                        mod = true;
                    case x12saDir{2}
                        fpath = fullfile(fpath, x12saDir{2});
                        mod = true;
                end
            end
            
            % no modifications to the path -> file not found
            if ~mod
                fpath = [];
                return
            end
                
        else
            if numel(nexusFiles)>1
                warning('More than one file found for scan S%05d! Using the first entry.', scanNr)
            end
            fpath = fullfile(nexusFiles(1).folder, nexusFiles(1).name);
        end
    end
end
