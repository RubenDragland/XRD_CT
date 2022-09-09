%NEXUS_READ read nexus file
% 
% ** fpath              file path or path to a parent directory of the
%                       nexus directory
%
% *optional*
% ** scanNr             scan number
% ** orientation        overwrite nexus orientation; should be given as [transpose rot90]
% ** filter             request a specific entry of the nexus file; default: /
% ** range              optional range to read only a portion of a datasets (cf. io.HDF.hdf5_load)
% ** attr               load HDF5 attributes; default: false
% ** mask               file path of a valid pixel mask; either .h5 or .mat
%
% returns:
% ++ nxs                requested nexus data, either as struct or cell of
%                       structs
% ++ mask               mask; only available if a mask file was specified
%
% EXAMPLES:
%       % specify the full path
%       nxs = io.nexus_read('~/Data10/data/S00000-00999/S00250/e16812_1_00250.h5');  
%
%       % base path + scanNr; this will load everything
%       nxs = io.nexus_read('~/Data10', 'scanNr', 250);
%
%       % load data entry ('/entry/data')
%       nxs = io.nexus_read('~/Data10', 'scanNr', 250, 'filter', 'data');
%
%       % load spec data
%       nxs = io.nexus_read('~/Data10', 'scanNr', 250, 'filter', 'spec');
%
%       % load pilatus_1 data
%       nxs = io.nexus_read('~/Data10', 'scanNr', 250, 'filter', 'pilatus_1');
%
%       % load pilatus_1 and pilatus_2
%       nxs = io.nexus_read('~/Data10', 'scanNr', 250, 'filter', 'pilatus');
%
%       % load all detectors
%       nxs = io.nexus_read('~/Data10', 'scanNr', 250, 'filter', 'detectors');
%
%       % load mcs data
%       nxs = io.nexus_read('~/Data10', 'scanNr', 250, 'filter', 'mcs');
%
%       % load sgalil data
%       nxs = io.nexus_read('~/Data10', 'scanNr', 250, 'filter', 'sgalil');
%
%       % load pilatus_1 + mask
%       [nxs, mask] = io.nexus_read('~/Data10', 'scanNr', 250, 'filter', 'pilatus_1', 'mask', '~/Data10/fancyMask.h5');
%
%       % load mask from an h5 mask file; use 'orientation' if you would like to adjust the orientation of the loaded mask
%       mask = io.nexus_read('./myFancyMask.h5', 'filter', 'mask');
%
%       % optional filters are inter alia: 
%       % spec, detectors, meta, energy, lambda, date, scan, eiger, pilatus, falcon...
%
%
% see also: io.HDF.hdf5_load()

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

function [nxs, maskOut] = nexus_read(fpath, varargin)

par = inputParser;
par.addParameter('scanNr', -1, @isnumeric)  
par.addParameter('orientation', [], @isnumeric)
par.addParameter('filter', [], @ischar)
par.addParameter('range', [], @iscell)
par.addParameter('attr', false, @islogical)
par.addParameter('mask', [], @ischar)

par.parse(varargin{:})
vars = par.Results;

vars.loc = [];
nxs = [];
maskOut = [];

% using '-a' and a range is currently not supported by hdf5_load
if ~isempty(vars.range) && vars.attr
    error('Loading attributes and reading a slice is currently not supported.')
end

% prepare hdf5_load options (cf. io.HDF.hdf5_load)
if ~isempty(vars.range)
    vars.hdf5_load_options = vars.range;
elseif vars.attr
    vars.hdf5_load_options = '-a';
else
    vars.hdf5_load_options = [];
end

% find nexus file
fpath = utils.find_nexus_file(fpath, vars.scanNr);
if isempty(fpath)
    error('Could not find nexus file.')
end
vars.fpath = fpath;


if isempty(vars.filter)
    vars.loc = '/';
    vars.applyOrientation = true;
    
elseif strcmpi(vars.filter(1), '/')
    vars.loc = vars.filter;
else
    % check version tag and beamline
    switch lower(vars.filter)  
        case {'spec', 'speces1'}
            vars.loc = '/entry/collection/specES1';
        case {'data'}
            vars.loc = '/entry/data';
            vars.applyOrientation = true;
        case {'detector', 'detectors'}
            % find all detectors
            vars.applyOrientation = true;
            h = h5info(fpath, '/entry/instrument');
            for ii=1:length(h.Groups)
                if isfield(h.Groups(ii).Attributes, 'Name') && strcmpi(h.Groups(ii).Attributes.Name, 'NX_class') && strcmpi(h.Groups(ii).Attributes.Value, 'NXdetector')
                    vars.loc{end+1} = h.Groups(ii).Name;
                end
            end
            
        case {'meta', 'collection'}
            vars.loc = '/entry/collection';
        case {'energy', 'mokev'}
            vars.loc = '/entry/collection/specES1/mokev';
        case {'lambda', 'wavelength'}
            vars.loc = '/entry/instrument/monochromator/wavelength';
            vars.scale = 1e-10;
        case {'bpm4', 'bpm'}
            vars.loc = '/entry/collection/specES1/bpmi';
        case {'time', 'd', 'date'}
            vars.loc = '/entry/collection/specES1/D';
        case {'scan', 'type', 's'}
            vars.loc = '/entry/collection/specES1/S';
        case 'eiger'
            % find all eigers
            vars.applyOrientation = true;
            h = h5info(fpath, '/entry/instrument');
            for ii=1:length(h.Groups)
                if isfield(h.Groups(ii).Attributes, 'Name') && strcmpi(h.Groups(ii).Attributes.Name, 'NX_class') && strcmpi(h.Groups(ii).Attributes.Value, 'NXdetector')
                    if contains(h.Groups(ii).Name, 'eiger')
                        vars.loc{end+1} = h.Groups(ii).Name;
                    end
                end
            end
        case {'eiger_1', 'eiger1'}
            vars.applyOrientation = true;
            vars.loc = '/entry/instrument/eiger_1';
        case {'eiger_4', 'eiger4'}
            vars.applyOrientation = true;
            vars.loc = '/entry/instrument/eiger_4';
        case 'pilatus'
            % find all pilati (?)
            vars.applyOrientation = true;
            h = h5info(fpath, '/entry/instrument');
            for ii=1:length(h.Groups)
                if isfield(h.Groups(ii).Attributes, 'Name') && strcmpi(h.Groups(ii).Attributes.Name, 'NX_class') && strcmpi(h.Groups(ii).Attributes.Value, 'NXdetector')
                    if contains(h.Groups(ii).Name, 'pilatus')
                        vars.loc{end+1} = h.Groups(ii).Name;
                    end
                end
            end
        case {'pilatus_1', 'pilatus1'}
            vars.applyOrientation = true;
            vars.loc = '/entry/instrument/pilatus_1';
        case {'pilatus_2', 'pilatus2'}
            vars.applyOrientation = true;
            vars.loc = '/entry/instrument/pilatus_2';
        case {'falconx1', 'falcon'}
            vars.applyOrientation = false;
            vars.loc = '/entry/instrument/FalconX1';
        case {'mcs', 'mcs_data'}
            vars.loc = '/entry/collection/mcs/mcs_data';
        case {'sgalil'}
            vars.loc = '/entry/collection/sgalil';
        otherwise
            % try to find the specified entry
            h = h5info(fpath, '/entry');
            [vars.loc, vars.applyOrientation] = find_loc(h, vars.filter);

    end
end


try
    % load file with specified options (vars)
    [nxs, vars] = read_loc(vars);
catch ME
    % return empty file if loading failed
    nxs = [];
end

if ~isempty(vars.mask) && ~isempty(nxs)
   [~, ~, ext] = fileparts(vars.mask);
   switch ext
       case '.mat'
           mask = load(vars.mask);
           maskOut = mask.mask;
           if isfield(nxs, 'data') && (size(maskOut,1)~= size(nxs.data,1) || size(maskOut,2)~=size(nxs.data,2))
               warning('Mask file dimensions and nexus data don''t match.')
           end
       case '.h5'
           maskOut = io.nexus_read(vars.mask, 'filter', 'mask', 'orientation', vars.orientation);
       otherwise
           error('Unknown file extension.')
   end
   
   % 
    
end


end

% read data from nexus file
function [nxs, vars] = read_loc(vars)

if iscell(vars.loc) && numel(vars.loc)>1
    for ii=1:numel(vars.loc)
        [nxs{ii}, vars] = load_data(vars.fpath, vars.loc{ii}, vars);
    end
elseif iscell(vars.loc) && numel(vars.loc)==1
    [nxs, vars] = load_data(vars.fpath, vars.loc{1}, vars);
else
    [nxs, vars] = load_data(vars.fpath, vars.loc, vars);
end

end

function [nxs, vars] = load_data(fpath, loc, vars)
    h = h5info(fpath, loc);
    isDetector = false;
    ii=1;
    while ii<=length(h.Attributes)
        if strcmpi(h.Attributes(ii).Name, 'NX_class') && strcmpi(h.Attributes(ii).Value, 'NXdetector')
            isDetector = true;
            break;
        end
        ii = ii+1;
    end
    
    % load data
    if isDetector
        
        arg.data_path = utils.abspath(fpath);
        arg.nthreads = 14;
        arg.precision = 'single';
        arg.extension = 'h5';
        arg.data_location = fullfile(loc, 'data');
        if isempty(vars.hdf5_load_options)
            try
                nxs.data = io.read_measurement(arg);
            catch
                warning('Fast loading routine failed.')
                nxs.data = io.HDF.hdf5_load(fpath, fullfile(loc, 'data'), vars.hdf5_load_options);
            end
        else
            nxs.data = io.HDF.hdf5_load(fpath, fullfile(loc, 'data'), vars.hdf5_load_options);
        end
    
        % load the rest
        for ii=1:length(h.Groups)
            locName = strsplit(h.Groups(ii).Name, '/');
            nxs.(locName{end}) = io.HDF.hdf5_load(fpath, h.Groups(ii).Name, vars.hdf5_load_options);
        end
        for ii=1:length(h.Datasets)
            if strcmpi(h.Datasets(ii).Name, 'data')
                continue
            end
            nxs.(h.Datasets(ii).Name) = io.HDF.hdf5_load(fpath, fullfile(loc, h.Datasets(ii).Name), vars.hdf5_load_options);
        end
        
        % load detector orientation
        if isempty(vars.orientation) && isfield(nxs, 'orientation')
            if vars.attr
                vars.orientation = [nxs.orientation.transpose.Value nxs.orientation.rot90.Value];
            else
                vars.orientation = [nxs.orientation.transpose nxs.orientation.rot90];
            end
        end

    else
        nxs = io.HDF.hdf5_load(fpath, loc, vars.hdf5_load_options);
    end
    
    if ~isempty(vars.orientation)
        if isDetector
            if vars.attr
                nxs.data.Value = math.applyTransform(nxs.data.Value,vars.orientation,'transpose_rot90');
            else
                nxs.data = math.applyTransform(nxs.data,vars.orientation,'transpose_rot90');
            end
        else
            if vars.attr
                nxs.Value = math.applyTransform(nxs.Value,vars.orientation,'transpose_rot90');
            else
                nxs = math.applyTransform(nxs,vars.orientation,'transpose_rot90');
            end
        end
    end
    
    if isfield(vars, 'scale')
        nxs = nxs.*vars.scale;
    end
    if isstruct(nxs)
        nxs.NXlocation = loc;
    end
    

end


% find target in nexus file 
function [loc, applyOrientation] = find_loc(h, locTarget)
    loc = [];
    applyOrientation = false;
    ii = 1;
    if isfield(h, 'Groups')
        while ii<=length(h.Groups)
            grps = strsplit(h.Groups(ii).Name, '/');
            if strcmpi(grps{end}, locTarget)
                loc = fullfile(h.Name, locTarget);
                if isfield(h.Groups(ii).Attributes, 'Name') && strcmpi(h.Groups(ii).Attributes.Name, 'NX_class') && strcmpi(h.Groups(ii).Attributes.Value, 'NXdetector')
                	applyOrientation = true;
                end
            else
                [loc, applyOrientation] = find_loc(h.Groups(ii), locTarget);
            end
            if ~isempty(loc)
                return;
            else
                ii = ii+1;
            end
        end
    end
    if isfield(h, 'Datasets')
        while ii<=length(h.Datasets)
            if strcmpi(h.Datasets(ii).Name, locTarget)
                loc = fullfile(h.Name, locTarget);
            else
                [loc, applyOrientation] = find_loc(h.Datasets(ii), locTarget);
            end
            if ~isempty(loc)
                return;
            else
                ii = ii+1;
            end
        end
    end
    if isfield(h, 'Attributes')
        while ii<=length(h.Attributes)
            if strcmpi(h.Attributes(ii).Name, locTarget)
                loc = fullfile(h.Name, locTarget);
            end
            if ~isempty(loc)
                return;
            else
                ii = ii+1;
            end
        end
    end
end


