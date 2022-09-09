%CONVERT2NEXUS converts and merges data into a nexus-formatted HDF5 file 
%
% *optional*
% ** scanNr             start at the given scan number
% ** finished           convert and ingest the last dataset, as well as the
%                       derived data
%
% EXAMPLES:
%   % start at scan number 1:
%       convert2NEXUS();
%
%   % start at scan number 150:
%       convert2NEXUS('scanNr', 150);
%
%   % wrap up the beamtime:
%       convert2NEXUS('finished', true);
%
% Pleas note that the script is designed to be used during an ongoing
% measurement, and therefore only converts n-1 datasets, that is it waits
% until the next measurement has started. Once your measurements are done,
% use convert2NEXUS('finished', true) to convert and ingest the last
% dataset.
% 
% If you want to rerun the converter on the ra cluster, you have to adjust 
% a few parameters. The following changes are needed to convert data 
% located in /das/work/p17/p171234/online_data:
%
%   % the eaccount cannot be detected automatically on ra. Please change it manually
%   s.eaccount = 'e171234';     
%
%   % the base path specifies the path to the software directories as well
%   as the path to the output directory. You might want to have a different 
%   output directory but make sure that the directory exists 
%   before you start the script!
%   s.base_path = '/das/work/p17/p171234/offline_data/';      
%
%   % if you are not working with the cxs_software file structure 
%   % (e.g. [s.base_path 'cxs_software/base']), set the
%   % beamline_file_structure to false. It will look for cSAXS_matlab_base
%   % and cSAXS_shell_scripts directories in s.base_path
%   s.beamline_file_structure = false;
%
%
%   % the raw_data_path defines the path to the detector directories
%   s.raw_data_path = '/das/work/p17/p171234/online_data/'
%
%   % disable online archiving
%   s.ingestData = false;
%
%   % make sure to double-check all "s.(detector-name).delete_raw_data" 
%   % flags and adjust them if needed.

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

function convert2NEXUS(varargin)
import utils.*
import beamline.nexus.*

check_input = @(x) islogical(x) || isnumeric(x);
par = inputParser;
par.addParameter('scanNr', 1, @isnumeric)               % starting point
par.addParameter('finished', false, check_input)        % !!!WIP!!! convert the last dataset and ingest everything

par.parse(varargin{:})
vars = par.Results;

% starting point
s.scanNr = vars.scanNr;

% general settings
s.verbose_level = 1;                                    % increase verbosity for more output; 1 is usually enough
s.debug = false;                                        % if debug==true, no system calls are performed
s.eaccount = beamline.identify_eaccount();              % when executed at the beamline, beamline.identify_eaccount() can be used instead of hardcoding e.g. 'e17970'
s.finished = vars.finished;                             % convert the last dataset and ingest everything
s.skipNonDatasets = true;                               % skip scans without dataset_id; this usually only happens for aborted scans. For legacy, to convert old datasets for which the dataset_id does not exist set to = true, note that the automatic archiving in petabyte archive will not work in either case

% input
s.beamline_file_structure = true;                       % file structure is the same as at the cSAXS beamline, that is the directories bin, cxs_software etc exist.
s.base_path = '~/Data10/';                               % e.g. ~/Data10/
s.raw_data_path = '~/Data10/';                           % e.g. ~/Data10/

% output
s.targetDir = 'data';                                   % output directory; could also be a detector directory
s.targetDirIsFlat = false;                              % set to false for creating a directory for each scan
s.fnamePrefix = [s.eaccount '_1'];                      % file name prefix; e.g. e18041_1
s.fnamePattern = sprintf('%s_%%05d.h5', s.fnamePrefix); % file name pattern; needs a placeholder for the scan number (e.g."%05d"); e.g. e18041_1_%05d.h5, which would result in e18041_1_00150.h5 for scan number 150.



%%%%%%%%%%%%%%%%%%%%%%
%%% raw data input %%%
%%%%%%%%%%%%%%%%%%%%%%
if s.beamline_file_structure
    bin_subdir = '/bin';
    base_subdir = 'cxs_software/base';
else
    bin_subdir = '/cSAXS_shell_scripts/beamline';
    base_subdir = 'cSAXS_matlab_base';
end


% Eiger1.5M
s.eiger1p5M.converter_path = fullfile(s.base_path, bin_subdir, 'eiger1p5M_converter/hdf5MakerOMNY');   % path to the converter binaries
s.eiger1p5M.path = fullfile(s.raw_data_path,'eiger_4');                             % path to the Eiger1.5M raw data
s.eiger1p5M.delete_raw_data = true;                                                 % delete raw data after conversion to HDF5

% Eiger500K
s.eiger500K.path = fullfile(s.raw_data_path,'eiger');                             % path to the Eiger500K raw data
s.eiger500K.delete_raw_data = false;                                                 % delete raw data after conversion

% Pilatus 2M
s.pilatus2m.path = fullfile(s.raw_data_path,'pilatus_1');                           % path the Pilatus raw data
s.pilatus2m.converter = fullfile(s.base_path, base_subdir, '/+beamline/+nexus/private/cbf2hdf5/cbf2hdf5'); % path to the cbf2hdf5 converter
s.pilatus2m.delete_raw_data = false;                                                % delete the raw data after conversion

% Pilatus 300K
s.pilatus300k.path = fullfile(s.raw_data_path,'pilatus_2');                         % path the Pilatus raw data
s.pilatus300k.converter = fullfile(s.base_path, base_subdir, '/+beamline/+nexus/private/cbf2hdf5/cbf2hdf5'); % path to the cbf2hdf5 converter
s.pilatus300k.delete_raw_data = false;                                              % delete the raw data after conversion

% Falcon
s.falcon.path = fullfile(s.raw_data_path,'falcon');                                 % path to the Falcon raw data
s.falcon.delete_raw_data = false;                                                   % delete the raw data after conversion

% SPEC data
s.spec.xmlLayoutFile = fullfile(s.base_path, bin_subdir, '/nexus/layout.xml');      % path to the nexus XML layout file
s.spec.specParser = fullfile(s.base_path, base_subdir, '/+io/spec_reader/spec_reader');  % path to the spec_reader binary
s.spec.appendEmpty = true;                                                          % write spec data, even if no other raw data exists

% orchestra
s.orchestra.path = fullfile(s.raw_data_path, '/specES1/scan_positions/');           % path to the orchestra files

% MCS
s.mcs.path = fullfile(s.raw_data_path, '/mcs');                                     % path to the mcs directory
s.mcs.entry = 'mcs_data';                                                           % dataset name

% Sgalil
s.sgalil.path = fullfile(s.raw_data_path, '/sgalil');                               % path to the sgalil directory

% archiving to CSCS
s.ingestion.cSAXSIngest = fullfile(s.base_path, bin_subdir, '/scicat/cSAXSIngest'); % path to cSAXSIngest binary
s.ingestion.tmpDir = '~/Documents/';                                                % output directory for intermediate files; ideally outside Data10
s.ingestion.ingestData = true;                                                     % !!!WIP!!! set to true for ingesting and archiving with CSCS;
s.ingestion.exclude = {'eiger_4', 'pilatus_1', 'pilatus_2'};                                     % exclude these directories

% queue the functions
s.funcSequence = {...
    @appendEiger1p5M, ...
    @appendEiger500K, ...
    @appendPilatus2M, ...
    @appendPilatus300K, ...
    @appendFalcon, ...
    @appendSpecData, ...
    @ingestData, ...
    @scanFinished, ...
    @deleteData ...
    };


% let's go!
beamline.nexus.startConversion(s);

end


