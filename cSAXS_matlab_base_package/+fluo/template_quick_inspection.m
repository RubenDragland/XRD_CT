% template_quick_inspection.m  
%                                                                                                   
% to process fluorescence data during the experiment

% This code and subroutines are part of a continuous development. There is no 
% liability on PSI or cSAXS. Before publishing results using this code 
% please notify M. Guizar-Sicairos, (mguizar@gmail.com)
% or another PSI cSAXS staff member.

close all
%% Prepare substructures for paths, processing and detector calibration 
%  this section is likely to be changed only once in the beginning of the 
%  experiment.

% paths substructure
fls.    path.base_path = '../../../online_data/';
fls.    path.matlab = '../../training_ana/';
fls.    path.xray_database = './+fluo/xraydb.mat'; % this is the database with X-ray emission and absorption lines
fls.    path.output = '../../analysis_fluo/';
fls.    path.marked_lines = []; % path to file with fluorescence lines that have been marked and identified on the spectrum
fls.    path.calibration = []; % file with detector calibration

% processing substructure
fls.    proc.metadata_source = 'orchestra'; % options are 'spec', 'orchestra' and 'sgalil'
fls.    proc.file_type = 'raw'; % options are 'nexus' and 'raw'
fls.    proc.motorx = []; % x-axis motor NOTE: ignored if the metadata_source = 'sgalil' or 'orchestra'
fls.    proc.motory = []; % y-axis motor NOTE: ignored if the metadata_source = 'sgalil' or 'orchestra'

fls.    proc.dead_time_correction = false; % as default better always correct the dead time
fls.    proc.stitched = 'none'; % options are 'none', 'blocks', 'lines'
fls.    proc.use_interpolation = false; % option for using interpolation when mapping the fluorescence data
fls.    proc.mapping_mode = true; % if mapping mode is true, then 2D map is assumed, otherwise just multiple spectra like in loopscan

% other parameters
fls.    calib.scan_num = [327]; % scans for calibration, can be a range of scans or a single scan
                                 % range of scans is summed up in a single spectrum, but if the scans represent
                                 % different projections instead of the stitches, only one scan will be
                                 % displayed during the calibration
fls.    calib.tested = false;

addpath(fls.path.matlab)
fls     = fluo.unpack_xray_database(fls);

%% Load or perform detector calibration using 2D map
%  This section allows you to calibrate the detector from channel numbers
%  into the units of keV. If the fls.calib.file is provided then you will
%  be invited to inspect whether the detector is properly calibrated. 
%  However, if there is no fls.calib.file provided then you will be asked 
%  to assist the script in calibrating the spectrum summed from the scans
%  given by fls.calib.scan_num.

fls = fluo.calibrate_detector(fls);
fls = fluo.test_detector_calibration(fls);

%% Inspect data
%  This section allows you to load specific scan and inspect specific
%  energy lines in this scan
fls.proc.mapping_mode = true; % if mapping mode is true, then 2D map is assumed, otherwise just multiple spectra like in loopscan
scans_to_inspect = 330;

fluo_data = fluo.fluo_read(fls.path.base_path, 'scanNr', scans_to_inspect, 'fluo_structure', fls);

if fls.proc.mapping_mode
    fluo.inspect_projection(fluo_data, fls)
else 
    fluo_data = fluo_data{1};
    if isfield(fls.calib,'energy')
        energy = fls.calib.energy;
        channels = false;
    else
        energy = [1:numel(spectrum)];
        channels = true;
    end
    spectrum = sum(fluo_data, 2);
    
    if channels
        semilogy(spectrum)
        xlabel('Channel number, [ ]')
    else
        semilogy(energy, spectrum)
        xlabel('Energy, [keV]')
    end
end


%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
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
%   using the “cSAXS matlab package” developed by the CXS group,
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