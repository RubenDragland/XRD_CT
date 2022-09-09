%FLUO_READ_STITCHED_EXTENDED read fluo file
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
% ++ real_time_data         map of real time values
% ++ live_time_data         map of live time values
% ++ correction_coeff_data  map of calculated correction coefficient
%
% see also: fluo.fluo_read()

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

function [fluo_data, real_time_data, live_time_data, correction_coeff_data] = fluo_read_stitched_extended(fpath, varargin)

par = inputParser;
par.addParameter('scanNr', -1, @isnumeric)
par.addParameter('file_type', [], @ischar)
par.addParameter('metadata_source', [], @ischar)
par.addParameter('dead_time_correction', true, @islogical)


par.parse(varargin{:})
vars = par.Results;

if isempty(vars.file_type)
    error('Please specify a valid file type.')
end

if isempty(vars.metadata_source)
    error('Please specify a valid metadata source.')
end

fluo_data = {};
real_time_data = {};
live_time_data = {};
correction_coeff_data = {};
x_loc = {};
y_loc = {};

for ii = 1 : numel(vars.scanNr)
    if strcmp(vars.file_type, 'nexus')
        fpath_loc = utils.find_nexus_file(fpath, vars.scanNr(ii));
        fprintf('Loading %s\n', fpath_loc)
        
        while ~utils.nexus_file_is_ready(fpath_loc)
            fprintf('%s is not ready. pause for 2 seconds\n', fpath_loc);
            pause(2);
        end
        
        pos_data_loc = io.nexus_read(fpath_loc, 'filter', vars.metadata_source);
        fluo_data_loc = io.nexus_read(fpath_loc, 'filter', 'falcon');
        if isempty(fluo_data_loc)
            error(sprintf('No falcon data found in the nexus file for scan %d.', vars.scanNr(ii)))
        end
        
    elseif strcmp(vars.file_type, 'raw')
        fpath_loc = find_raw_fluo_file(fpath, vars.scanNr(ii));
        if ~isempty(fpath_loc)
            fprintf('Loading %s\n', fpath_loc)
            scanNr = strsplit(fpath_loc, '/');
            scanNr = scanNr{end-1};
            scanNr = str2num(scanNr(2:end));
        
            fpath_pos_data_loc = find_pos_data_file(fpath_loc, scanNr);
            pos_data_loc = beamline.read_omny_pos(fpath_pos_data_loc);
            fluo_data_loc = io.HDF.hdf5_load(fpath_loc);
        end
    else
        error('Not implemented data type or wrong spelling of the data type.')
    end
    
    if strcmp(vars.metadata_source, 'orchestra')
        x_loc{ii} = double(-pos_data_loc.Average_x_st_fzp);
        y_loc{ii} = double(-pos_data_loc.Average_y_st_fzp);
    elseif strcmp(vars.metadata_source, 'spec')
        x_loc{ii} = double(-pos_data_loc.px);
        y_loc{ii} = double(-pos_data_loc.py);
    else
        error('Not implemented metadata source.');
    end
    
    vars.number_of_spectra = numel(x_loc{ii});
    if isempty(fpath_loc)
        fluo_data_loc = zeros(2048, 1, numel(x_loc{ii})) .* NaN; % HARD CODED NUMBER OF CHANNELS!!!
    else
        [fluo_data_loc, time] = fluo_read_loc(fluo_data_loc, vars.dead_time_correction, vars.number_of_spectra);
    end
    
    [~, num_of_detectors, ~] = size(fluo_data_loc);
    if num_of_detectors == 1
        fluo_data_loc = squeeze(fluo_data_loc);
    else
        error('Currently fluo data reader is only implemented for 1 detector, but %d are given', num_of_detectors)
    end
    
    fluo_data{ii} = fluo_data_loc;
    real_time_data{ii} = time.real_time;
    live_time_data{ii} = time.live_time;
    correction_coeff_data{ii} = time.correction_coef;
end


fluo_data = [fluo_data{:}];
real_time_data = [real_time_data{:}];
live_time_data = [live_time_data{:}];
correction_coeff_data = [correction_coeff_data{:}];
x_loc_loc = [];
y_loc_loc = [];
for kk = 1 : numel(x_loc)
    x_loc_loc = [x_loc_loc x_loc{kk}'];
    y_loc_loc = [y_loc_loc y_loc{kk}'];
end
x_loc = x_loc_loc;
y_loc = y_loc_loc;
%scatter(x_loc, y_loc)
[fluo_data, real_time_data, live_time_data, correction_coeff_data] = fluo_reshape_loc(fluo_data, real_time_data, live_time_data, correction_coeff_data, x_loc, y_loc);

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

function [fluo_data_loc, time] = fluo_read_loc(fluo_data_loc, dead_time_correction, number_of_spectra)
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


live_time = decoded_data.LiveTime(1:number_of_spectra);
real_time = decoded_data.RealTime(1:number_of_spectra);
correction_coef = real_time ./ live_time;
if dead_time_correction == 1
    for ch_num = 1:number_of_spectra
        fluo_data_loc(:, 1, ch_num) = fluo_data_loc(:, 1, ch_num) * correction_coef(ch_num);
    end
end
time.live_time = live_time;
time.real_time = real_time;
time.correction_coef = correction_coef;

end

function [fluo_data_loc, real_time_data, live_time_data, correction_coeff_data] = fluo_reshape_loc(fluo_data_loc, real_time_data, live_time_data, correction_coeff_data, x_loc, y_loc)
% threshold = mean(abs(diff(sort(x_loc)))) / 2;
% num_of_cols = abs(diff(sort(x_loc))) > threshold;
% Nx = sum(num_of_cols);
% 
% threshold = mean(abs(diff(sort(y_loc)))) / 2;
% num_of_rows = abs(diff(sort(y_loc))) > threshold;
% Ny = sum(num_of_rows);
% 
% if mod(numel(x_loc), Nx) == 0 && mod(numel(x_loc), Ny) ~= 0
%     Ny = numel(x_loc) / Nx;
% elseif mod(numel(x_loc), Ny) == 0 && mod(numel(x_loc), Nx) ~= 0
%     Nx = numel(x_loc) / Ny;
% elseif mod(numel(x_loc), Nx) ~= 0 && mod(numel(x_loc), Ny) ~= 0
%     warning('Inconsistency in number of points or the estimated number of rows and columns./nTrying different method.')
%     
% end

Nx = floor(sqrt((max(x_loc) - min(x_loc)) * numel(x_loc) / (max(y_loc) - min(y_loc))));
Ny = round(numel(x_loc) / Nx);

[num_of_channels, num_of_spectra] = size(fluo_data_loc);

xg = linspace(min(x_loc), max(x_loc), Nx);
yg = linspace(min(y_loc), max(y_loc), Ny);
[Xg, Yg] = meshgrid(xg, yg);

intermid_data = zeros(num_of_channels, Ny, Nx);
for ee = 1:num_of_channels
    F = scatteredInterpolant(x_loc', y_loc', fluo_data_loc(ee, :).');
    intermid_data(ee, :, :) = F(Xg, Yg);
end
fluo_data_loc = permute(intermid_data,[1 2 3]);

F = scatteredInterpolant(x_loc', y_loc', real_time_data.');
intermid_data = F(Xg, Yg);
real_time_data = permute(intermid_data,[1 2]);

F = scatteredInterpolant(x_loc', y_loc', live_time_data.');
intermid_data = F(Xg, Yg);
live_time_data = permute(intermid_data,[1 2]);

F = scatteredInterpolant(x_loc', y_loc', correction_coeff_data.');
intermid_data = F(Xg, Yg);
correction_coeff_data = permute(intermid_data,[1 2]);

end