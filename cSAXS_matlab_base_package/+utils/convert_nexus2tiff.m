%CONVERT_NEXUS2TIFF extract data from a nexus file and save it as a tiff
%stack or mutliple tiff files
% 
% ** scanNr             scan number
% ** filter             nexus filter; e.g. 'pilatus_1'
%
% *optional*
% ** basePath           base path; default: '~/Data10/'
% ** outdir             output directory path; default: '~/Data10/tiffs'
% ** ext                target file extension; default: 'tif'
% ** writeStack         if true, a single tif file contains a 3D stack. Set
%                       to false for writing mutliple files. Default: true
% ** flatOutput         write into a single output directory, grouped into 
%                       packages of 1000, instead of a directory for each 
%                       file; default: false 
%
% see also: io.nexus_read

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


function convert_nexus2tiff(scanNr, filter, varargin)

par = inputParser;
par.addParameter('basePath', '~/Data10', @ischar)
par.addParameter('outdir', '~/Data10/tiffs/', @ischar)
par.addParameter('ext', 'tif', @ischar)
par.addParameter('writeStack', true, @islogical)
par.addParameter('flatOutput', false, @islogical)


par.parse(varargin{:})
vars = par.Results;

% get nexus filename
fname = utils.find_nexus_file(utils.abspath(vars.basePath), scanNr);

% prepare output directories
vars.outdir = fullfile(utils.abspath(vars.outdir), utils.compile_x12sa_dirname(scanNr, vars.flatOutput));
if ~exist(vars.outdir, 'dir')
    mkdir(vars.outdir);
end


% read nexus data and convert it to uint16
h = io.nexus_read(fname, 'filter', filter, 'orientation', [0 0]);
if isempty(h)
    error('Failed to read nexus file.');
end
h.data = uint16(h.data);

[~, fname_sub] = fileparts(fname);

if vars.writeStack
    % write tiffs in a single file
    fname_out = fullfile(vars.outdir, [fname_sub '.' vars.ext]);
    if exist(fname_out, 'file')
        delete(fname_out);
    end
    for ii=1:size(h.data,3)
        imwrite(h.data(:,:,ii), fname_out, 'tiff', 'WriteMode', 'append');
    end
else
    % write multiple files
    for ii=1:size(h.data,3)
        imwrite(h.data(:,:,ii), fullfile(vars.outdir, sprintf([fname_sub '_%05u.' vars.ext], ii)), 'tiff', 'WriteMode', 'overwrite');
    end
end

end