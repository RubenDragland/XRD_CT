%GET_DETECTOR_ORIENTATION read detector orientation from nexus file
% 
% ** detector       detector name, e.g. 'pilatus_1', 'eiger_4'...
% ** filename       full path of the nexus file from which you want to read
%
% *optional*
% ** orientType     type of the orientation vector. Either
%                   'transpose_rot90', 'flipud_rot90' or 
%                   'transpose_fliplr_flipud'; default: 'transpose_rot90'
%
% returns:
% ++ orientation    orientation vector
%
% EXAMPLES:
%       orientation = utils.get_detector_orientation('pilatus_1','./e17970_1_00494.h5');
%       orientation = utils.get_detector_orientation('pilatus_1','./e17970_1_00494.h5', 'orientType', 'transpose_fliplr_flipud');
% 
%       % loading cSAXS default orientation (not recommended)
%       orientation = utils.get_detector_orientation('pilatus_1');
%
% see also: math.applyTransform()

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



function orientation = get_detector_orientation(detector, varargin)

% process varargin
par = inputParser;
par.addParameter('orientType', [], @ischar)  
par.addParameter('filename', [], @ischar)  

par.parse(varargin{:})
vars = par.Results;

if nargin==0
    error('Please specify a detector.')
end

% get default path 
if strcmpi(detector(1),'/')
    detectorFullPath = detector;
else
    detectorFullPath = fullfile('/entry/instrument/', detector);
end

% try to read orientation from file
if ~isempty(vars.filename)
    % read orientation from nexus file
    d = [];
    try
        d = io.HDF.hdf5_load(vars.filename, fullfile(detectorFullPath, 'orientation'));
    catch
        if ~isempty(d) && isstruct(d)
            orientation = [d.transpose d.rot90];
        else
            warning('Failed to read orientation from file.')
            orientation = getcSAXSDefaults(detector);
            warning('Using cSAXS default orientation. Use with caution.')
        end
    end
    
else
    orientation = getcSAXSDefaults(detector);
    warning('Using cSAXS default orientation. Use with caution.')
end

if ~isempty(vars.orientType)
    orientation = math.convert_transform(orientation, 'transpose_rot90', vars.orientType);
end

end


function orientation = getcSAXSDefaults(detector)

if contains(detector, 'pilatus_1')
    orientation = [1 2];
elseif contains(detector, 'pilatus_2')
    orientation = [1 2];
elseif contains(detector, 'eiger_1')
    orientation = [1 0];
elseif contains(detector, 'eiger_4')
    orientation = [1 0];
elseif contains(detector, 'moench')
    orientation = [1 2];
else
    error('Unknown detector.')
end

end
    
        