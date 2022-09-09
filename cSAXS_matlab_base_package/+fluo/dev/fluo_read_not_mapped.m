%FLUO_READ read fluo file
%
% ** fpath                  file path or path to a parent directory
%
% *optional*
% ** scanNr                 scan number
% ** file_type              eiher 'raw' or 'nexus'
% ** metadata_source        source of positions, options are 'spec' and 'orchestra'
% ** dead_time_correction   option (true or false) for not correcting for the dead time
%                           the default is to correct (true)
%
% returns:
% ++ fluo_data              fluorescence data not mapped to 2D, here it is
%                           given in the shape of (points x channels)
%
% EXAMPLES:
%       % specify the full path and nexus as data type
%       fluo_data = fluo_read('~/Data10/data/S00000-00999/S00250/e16812_1_00250.h5', 'file_type', 'nexus', 'metadata_source', 'spec')
%
%       % base path + scanNr as range and file type
%       fluo_data = fluo_read('~/Data10', 'scanNr', 250, 'file_type', 'nexus');
%       % or
%       fluo_data = fluo_read('~/Data10', 'scanNr', 250, 'file_type', 'raw', 'metadata_source', 'spec');
%
% see also: fluo.falcon_decode_buffers()

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

function fluo_data = fluo_read_not_mapped(fpath, varargin)

par = inputParser;
par.addParameter('fluo_structure', [], @isstruct)
par.addParameter('scanNr', -1, @isnumeric)
par.addParameter('file_type', [], @ischar)
par.addParameter('metadata_source', [], @ischar)

par.parse(varargin{:})
vars = par.Results;

if ~isempty(vars.fluo_structure)
    file_type = vars.fluo_structure.proc.file_type;
    metadata_source = vars.fluo_structure.proc.metadata_source;
else
    file_type = vars.file_type;
    metadata_source = vars.metadata_source;
end

if isempty(file_type)
    error('Please specify a valid file type.')
end

if strcmp(file_type, 'nexus')
    fpath_loc = utils.find_nexus_file(fpath, vars.scanNr);
    fprintf('Loading %s\n', fpath_loc)

    while ~utils.nexus_file_is_ready(fpath_loc)
        fprintf('%s is not ready. pause for 2 seconds\n', fpath_loc);
        pause(2);
    end

    pos_data_loc = io.nexus_read(fpath_loc, 'filter', metadata_source);
    fluo_data_loc = io.nexus_read(fpath_loc, 'filter', 'falcon');
    if isempty(fluo_data_loc)
        error(sprintf('No falcon data found in the nexus file for scan %d.', vars.scanNr(ii)))
    end

elseif strcmp(file_type, 'raw')
    fpath_loc = find_raw_fluo_file(fpath, vars.scanNr);
    fprintf('Loading %s\n', fpath_loc)
    scanNr = strsplit(fpath_loc, '/');
    scanNr = scanNr{end-1};
    scanNr = str2num(scanNr(2:end));

    fpath_pos_data_loc = find_pos_data_file(fpath_loc, scanNr);
    pos_data_loc = beamline.read_omny_pos(fpath_pos_data_loc);
    fluo_data_loc = io.HDF.hdf5_load(fpath_loc);
else
    error('Not implemented data type or wrong spelling of the data type.')
end

if strcmp(metadata_source, 'orchestra')
    x_loc = double(-pos_data_loc.Average_x_st_fzp);
    y_loc = double(-pos_data_loc.Average_y_st_fzp);
elseif strcmp(metadata_source, 'spec')
    x_loc = double(-pos_data_loc.diode);
    y_loc = double(-pos_data_loc.diode);
else
    error('Not implemented metadata source.');
end

number_of_spectra = numel(x_loc);
if isempty(fpath_loc)
    fluo_data_loc = zeros(2048, 1, number_of_spectra) .* NaN; % HARD CODED NUMBER OF CHANNELS!!!
else
    fluo_data_loc = fluo_read_loc(fluo_data_loc, dead_time_correction, number_of_spectra);
end

[~, num_of_detectors, ~] = size(fluo_data_loc);
if num_of_detectors == 1
    fluo_data_loc = squeeze(fluo_data_loc);
else
    error('Currently fluo data reader is only implemented for 1 detector, but %d are given', num_of_detectors)
end

fluo_data{ii} = fluo_data_loc;
end

function [fpath] = find_raw_fluo_file(fpath, varargin)

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
    fluo_file = dir(fullfile(fpath, 'falcon.h5'));
    if isempty(fluo_file)
        mod = false; % keep track of modifications
        res = dir(fpath);
        for ii=1:length(res)
            switch res(ii).name
                case 'Data10'
                    fpath = fullfile(fpath, 'Data10');
                    mod = true;
                case 'falcon'
                    fpath = fullfile(fpath, 'falcon');
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
        if numel(fluo_file) > 1
            warning('More than one file found for scan S%05d! Using the first entry.', scanNr)
        end
        fpath = fullfile(fluo_file(1).folder, fluo_file(1).name);
    end
end

end

function [fpath] = find_pos_data_file(fpath, scanNr)
x12saDir = strsplit(fpath, '/');
fpath = [];
for path_iter = 1:numel(x12saDir)-4
    fpath = [fpath x12saDir{path_iter} '/'];
end

while exist(fpath, 'file')~=2
    % check current directory
    fluo_file = dir(fullfile(fpath, sprintf('*/scan_positions/*%u*.dat', scanNr)));
    if isempty(fluo_file)
        mod = false; % keep track of modifications
        res = dir(fpath);
        for ii=1:length(res)
            switch res(ii).name
                case 'specES1'
                    fpath = fullfile(fpath, 'specES1');
                    mod = true;
            end
        end
        
        % no modifications to the path -> file not found
        if ~mod
            fpath = [];
            return
        end
        
    else
        if numel(fluo_file) > 1
            warning('More than one file found for scan S%05d! Using the first entry.', scanNr)
        end
        fpath = fullfile(fluo_file(1).folder, fluo_file(1).name);
    end
end

end

function fluo_data_loc = fluo_read_loc(fluo_data_loc, dead_time_correction, number_of_spectra)
if isstruct(fluo_data_loc)
    try
        try
            buffer_data = fluo_data_loc.entry.data.data;
        catch ME
        end
        try
            buffer_data = fluo_data_loc.data;
        catch
        end
    catch ME
        error('Not implemented fluorescence structure. Falcon decode buffer accepts /entry/data/data or already extracted array.')
    end
else
    buffer_data = fluo_data_loc;
end

decoded_data = fluo.falcon_decode_buffers(buffer_data);
fluo_data_loc= decoded_data.Data(:, :, 1:number_of_spectra);

if dead_time_correction == 1
    live_time = decoded_data.LiveTime(1:number_of_spectra);
    real_time = decoded_data.RealTime(1:number_of_spectra);
    correction_coef = real_time ./ live_time;
    
    for ch_num = 1:number_of_spectra
        fluo_data_loc(:, ch_num) = fluo_data_loc(:, ch_num) * correction_coef(ch_num);
    end
end
end