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
% ++ fluo_data              requested fluorescence data
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

function fluo_data = fluo_read(fpath, varargin)

par = inputParser;
par.addParameter('fluo_structure', [], @isstruct)
par.addParameter('scanNr', -1, @isnumeric)
par.addParameter('file_type', [], @ischar)
par.addParameter('metadata_source', [], @ischar)
par.addParameter('dead_time_correction', true, @islogical)
par.addParameter('stitched', 'none', @ischar)
par.addParameter('use_interpolation', false, @islogical)
par.addParameter('motorx', [], @ischar)
par.addParameter('motory', [], @ischar)
par.addParameter('mapping_mode', true, @islogical)


par.parse(varargin{:})
vars = par.Results;

if ~isempty(vars.fluo_structure)
    file_type = vars.fluo_structure.proc.file_type;
    metadata_source = vars.fluo_structure.proc.metadata_source;
    dead_time_correction = vars.fluo_structure.proc.dead_time_correction;
    stitched = vars.fluo_structure.proc.stitched;
    use_interpolation = vars.fluo_structure.proc.use_interpolation;
    motorx = vars.fluo_structure.proc.motorx;
    motory = vars.fluo_structure.proc.motory;
    mapping_mode = vars.fluo_structure.proc.mapping_mode;
else
    file_type = vars.file_type;
    metadata_source = vars.metadata_source;
    dead_time_correction = vars.dead_time_correction;
    stitched = vars.stitched;
    use_interpolation = vars.use_interpolation;
    motorx = vars.motorx;
    motory = vars.motory;
    mapping_mode = vars.mapping_mode;
end

if isempty(file_type)
    error('Please specify a valid file type.')
end

if isempty(metadata_source)
    error('Please specify a valid metadata source.')
end

fluo_data = {};
x_loc = {};
y_loc = {};

f = waitbar(0, 'Loading fluorescence projections...');
for ii = 1 : numel(vars.scanNr)
    if strcmp(file_type, 'nexus')
        fpath_loc = utils.find_nexus_file(fpath, vars.scanNr(ii));
        check_nexus_file = utils.nexus_file_is_ready(fpath_loc);
        while ~check_nexus_file
            message_out = fprintf('Scan %d is not ready. pause for 2 seconds\n', vars.scanNr(ii));
            waitbar(ii/numel(vars.scanNr), f, message_out);
            pause(2);
            fpath_loc = utils.find_nexus_file(fpath, vars.scanNr(ii));
            check_nexus_file = utils.nexus_file_is_ready(fpath_loc);
        end
        message_out = fprintf('Loading scan number %d\n', vars.scanNr(ii));
        waitbar(ii/numel(vars.scanNr), f, message_out);
        
        pos_data_loc = io.nexus_read(fpath_loc, 'filter', metadata_source);
        fluo_data_loc = io.nexus_read(fpath_loc, 'filter', 'falcon');
        number_of_points = numel(io.nexus_read(fpath_loc, 'filter', 'Epoch'));
        number_of_bursts = io.nexus_read(fpath_loc, 'filter', 'burstn');
        number_of_spectra = number_of_points * number_of_bursts;
        if isempty(fluo_data_loc)
            message_out = sprintf('No falcon data found in the nexus file for scan %d.', vars.scanNr(ii));
            waitbar(ii/numel(vars.scanNr), f, message_out);
            error(message_out)
        end
        
    elseif strcmp(file_type, 'raw')
        fpath_loc = find_raw_fluo_file(fpath, vars.scanNr(ii));
        message_out = sprintf('Loading scan number %d', vars.scanNr(ii));
        waitbar(ii/numel(vars.scanNr), f, message_out);
        scanNr = strsplit(fpath_loc, '/');
        scanNr = scanNr{end-1};
        scanNr = str2num(scanNr(2:end));
        
        % !!!!!!!!!!!!!!!!!!!!! make a check of the next scan started for
        % the raw file case !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        path_to_spec = find_spec_file(fpath_loc);
        number_of_points = numel(io.spec_read(path_to_spec, 'ScanNr', scanNr).Epoch);
        number_of_bursts = io.spec_read(path_to_spec, 'ScanNr', scanNr).burstn;
        number_of_spectra = number_of_points * number_of_bursts;
        
        if strcmp(metadata_source, 'sgalil') || strcmp(metadata_source, 'orchestra')
            fpath_pos_data_loc = find_pos_data_file(fpath_loc, scanNr, metadata_source);
            pos_data_loc = beamline.read_omny_pos(fpath_pos_data_loc);
        elseif strcmp(metadata_source, 'spec')
            pos_data_loc = io.spec_read(path_to_spec, 'ScanNr', scanNr);
        else
            message_out = 'Not implemented metadata source';
            waitbar(ii/numel(vars.scanNr), f, message_out);
            error(message_out)
        end
        
        fluo_data_loc = io.nexus_read(fpath_loc,'filter','falcon');

        
    else
        message_out = 'Not implemented data type or wrong spelling of the data type.';
        waitbar(ii/numel(vars.scanNr), f, message_out);
        error(message_out)
    end
    
    if strcmp(metadata_source, 'orchestra')
        x_loc{ii} = double(-pos_data_loc.Average_x_st_fzp);
        y_loc{ii} = double(-pos_data_loc.Average_y_st_fzp);
    elseif strcmp(metadata_source, 'spec')
        x_loc{ii} = double(-pos_data_loc.(motorx));
        y_loc{ii} = double(-pos_data_loc.(motory));
    elseif strcmp(metadata_source, 'sgalil')
        x_loc{ii} = double(-pos_data_loc.Avg_x);
        y_loc{ii} = double(-pos_data_loc.Avg_y);
    else
        message_out = 'Not implemented metadata source';
        waitbar(ii/numel(vars.scanNr), f, message_out);
        error(message_out)
    end
    
    fluo_data_singlescan = fluo.fluo_raw_to_spectra(fluo_data_loc, 'dead_time_correction', dead_time_correction);
    clear fluo_read_loc
    
    
    [~, num_of_detectors, ~] = size(fluo_data_singlescan);
    if num_of_detectors == 1
        fluo_data_singlescan = squeeze(fluo_data_singlescan);
    else
        message_out = sprintf('Currently fluo data reader is only implemented for 1 detector, but %d are given', num_of_detectors);
        waitbar(ii/numel(vars.scanNr), f, message_out);
        error(message_out)
    end
    
    fluo_data{ii} = fluo_data_singlescan;  % MGS. This dynamic allocation of memory could slow down the reading quite a lot
end



if strcmp(stitched, 'none')
    for ii = 1 : numel(vars.scanNr)
        fluo_data{ii} = fluo_reshape_loc(fluo_data{ii}, x_loc{ii}, y_loc{ii}, stitched, use_interpolation, mapping_mode);
    end
elseif strcmp(stitched, 'block') || strcmp(stitched, 'lines')
    fluo_data = [fluo_data{:}];
    x_loc_loc = [];
    y_loc_loc = [];
    for kk = 1 : numel(x_loc)
        x_loc_loc = [x_loc_loc x_loc{kk}'];
        y_loc_loc = [y_loc_loc y_loc{kk}'];
    end
    x_loc = x_loc_loc;
    y_loc = y_loc_loc;
    fluo_data = fluo_reshape_loc(fluo_data, x_loc, y_loc, stitched, use_interpolation, mapping_mode);
else
    message_out = 'Not implemented stitching type';
    waitbar(ii/numel(vars.scanNr), f, message_out);
    error(message_out)
end

close(f)
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
        message_out = sprintf('File %s does not exist', fpath);
        waitbar(ii/numel(vars.scanNr), f, message_out);
        error(message_out)
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


function [fpath] = find_spec_file(fpath)


x12saDir = strsplit(fpath, '/');
fpath = [];
for path_iter = 1:numel(x12saDir)-4
    fpath = [fpath x12saDir{path_iter} '/'];
end

fpath = beamline.find_specDatFile(fpath);

end



function [fpath] = find_pos_data_file(fpath, scanNr, metadata_source)
x12saDir = strsplit(fpath, '/');
fpath = [];
for path_iter = 1:numel(x12saDir)-4
    fpath = [fpath x12saDir{path_iter} '/'];
end

if strcmp(metadata_source, 'orchestra')
    pos_folder_name = 'scan_positions';
elseif strcmp(metadata_source, 'sgalil')
    pos_folder_name = 'sgalil';
else
    error('Not implemented metadata source')
end

while exist(fpath, 'file')~=2
    % check current directory
    pos_file = dir(fullfile(fpath, sprintf('*/*%u*.dat', scanNr)));
    if isempty(pos_file)
        mod = false; % keep track of modifications
        res = dir(fpath);
        for ii=1:length(res)
            switch res(ii).name
                case 'specES1'
                    fpath = fullfile(fpath, 'specES1');
                    mod = true;
                case pos_folder_name
                    fpath = fullfile(fpath, pos_folder_name);
                    mod = true;
            end
        end
        
        % no modifications to the path -> file not found
        if ~mod
            fpath = [];
            return
        end
        
    else
        if numel(pos_file) > 1
            warning('More than one file found for scan S%05d! Using the first entry.', scanNr)
        end
        fpath = fullfile(pos_file(1).folder, pos_file(1).name);
    end
end

end

function fluo_data_loc = fluo_reshape_loc(fluo_data_loc, x_loc, y_loc, stitched, use_interpolation, mapping_mode)

if mapping_mode
    if strcmp(stitched, 'none') || strcmp(stitched, 'block')
        Nx = floor(sqrt((max(x_loc) - min(x_loc)) * numel(x_loc) / (max(y_loc) - min(y_loc))));
        Ny = round(numel(x_loc) / Nx);
    elseif strcmp(stitched, 'lines')
        error('Lines stitching is not yet implemented')
    else
        error('Not implemented stitching type')
    end

    [num_of_channels, ~] = size(fluo_data_loc);

    if ~use_interpolation
        fluo_data_loc = reshape(fluo_data_loc, [num_of_channels, Nx, Ny]);
        fluo_data_loc = permute(fluo_data_loc,[1 3 2]);
    else
        xg = linspace(min(x_loc), max(x_loc), Nx);
        yg = linspace(min(y_loc), max(y_loc), Ny);
        [Xg, Yg] = meshgrid(xg, yg);

        intermid_data = zeros(num_of_channels, Ny, Nx);
        for ii = 1:num_of_channels
            F = scatteredInterpolant(x_loc, y_loc, fluo_data_loc(ii, :).');
            intermid_data(ii, :, :) = F(Xg, Yg);
        end
        fluo_data_loc = intermid_data;
    end
end

end