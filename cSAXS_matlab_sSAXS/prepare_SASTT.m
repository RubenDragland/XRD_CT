% [projection] = prepare_SASTT(four_dat, data_process)
% function used for online loading of the data for SASTT
% fnames: structure with the names of the paths used for saving the data
% par: parameters for plotting and scan (par.snake, par.fast_axis)
% data_process: symmetric intensity data averaged over the segments for the
% whole q range
% four_dat: result of fft from analyze_one_scan
% output:
% projection: structure used as input for SASTT code

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
%
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing la this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS scanning SAXS package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% Additionally, any publication using the package, or any translation of the 
%     code into another computing language should cite:
%    O. Bunk, M. Bech, T. H. Jensen, R. Feidenhans'l, T. Binderup, A. Menzel 
%    and F Pfeiffer, “Multimodal x-ray scatter imaging,” New J. Phys. 11,
%    123016 (2009). (doi: 10.1088/1367-2630/11/12/123016)
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [projection] = prepare_SASTT(four_dat, data_process)

fnames = four_dat.fnames;
par = four_dat.par;


%determine the directory and file name
dirname = fullfile(fnames.save_dir, 'SASTT', fnames.sub_dir, 'projection_data');
filename = fullfile(dirname, sprintf('SASTT_%s.mat',par.sample_name));

%make a directory if does not exist
if ~exist(dirname, 'dir')
    mkdir(dirname);
end

% start an empty structure
projection = struct();
% organize the data into the structure
projection.first_scan_rot = (fnames(1).first_scan_no);
if (~isempty(fnames.basedir_trans_data))
    projection.diode = single(four_dat(1).trans_data);
end
%reshape the data matrix
if isfield(data_process,'sector_data')
    test1 = data_process.sector_data;
else
    test1 = [];
end
if isfield(data_process,'sector_data_q')
    data_sector_qres = data_process.sector_data_q;
else
    data_sector_qres=[];
end

%%% Adjust projection based on snake_scan or fast_axis %%%
if isfield(four_dat,'original_pos_data')
    pos_data = four_dat.original_pos_data;
else
    pos_data = [];
end
test1            = utils.adjust_projection(test1, par.snake_scan, par.fast_axis_x, pos_data);
data_sector_qres = utils.adjust_projection(data_sector_qres, par.snake_scan, par.fast_axis_x, pos_data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isempty(fnames.basedir_trans_data))
    %save and normalize by the transmission
    projection.data = single(test1./four_dat(1).trans_data);
    if ~isempty(data_sector_qres)
        projection.data_q = single(data_sector_qres./four_dat(1).trans_data);
    end
else
    projection.data = single(test1);
    if ~isempty(data_sector_qres)
        projection.data_q = single(data_sector_qres);
    end
    warning('There is no normalization with the transmission for SASTT. Beware.')
end

%read angle values
if par.use_nexus
    spec_data = io.nexus_read(fnames.base_dir_beamtime,'ScanNr',fnames(1).first_scan_no,'filter',fnames.spec_file);
else
    spec_data = io.spec_read(fnames.spec_file, 'ScanNr', fnames(1).first_scan_no);
end
rot_x = spec_data.rotx;
rot_y = spec_data.roty;
projection.rot_x = (rot_x);
projection.rot_y = (rot_y);

%calculate rotation matrix
if par.tomo_axis_x
    %Rotation of object around y: tilt axis B
    rotangle = deg2rad(rot_y); % transform to radians
    R = [cos(rotangle)  0   sin(rotangle);
        0         1         0;
        -sin(rotangle)  0   cos(rotangle)];
    %Rotation of object around x: tomo axis A
    rotangle = deg2rad(rot_x); % transform to radians
    R = R * [1         0             0;
        0   cos(rotangle)  -sin(rotangle);
        0   sin(rotangle)  cos(rotangle)];
    %save for rotation matrix for next steps
    projection.Rot_exp = single(R);
else % if tomo_axis_x = 0
    %Rotation of object around x: tilt axis A
    rotangle = deg2rad(rot_x); % transform to radians
    R =     [1      0               0;
        0  cos(rotangle)   -sin(rotangle);
        0  sin(rotangle)   cos(rotangle)];
    %Rotation of object around y: tomo axis B
    rotangle = deg2rad(rot_y); % transform to radians
    R =R * [cos(rotangle) 0  sin(rotangle);
        0      1        0      ;
        -sin(rotangle) 0  cos(rotangle)];
    
    %save for rotation matrix for next steps
    projection.Rot_exp = (R);
end
projection.par           = par;
projection.fnames        = fnames;
projection.integ         = four_dat.integ;
projection.scan_num      = four_dat.scan_num;
projection.scan_point    = four_dat.scan_point;
projection.positions_out = four_dat.positions_out;

% rename the projection structure to save as different variables
save_name = sprintf('p%d', (fnames(1).first_scan_no));
proj.(save_name) = projection;

% save or append to the SASTT input file
% this is used to avoid loading and saving while measuring/processing the
% input for SASTT. This works for Matlab v2017b or higher version

if ~exist(filename, 'file') % save first projection
    save(filename, '-struct', 'proj')
    display(['Saving SASTT projection to: ' filename])
else % append to the already existent file
    save(filename, '-struct', 'proj', '-append')
    display(['Appending SASTT projection to: ' filename])
end


