%INTEGRATE_RANGE_NEXUS  Azimuthal integration of SAXS patterns from NeXus
%   or H5 files
%
% beamline.integrate_range_nexus(startScanNr, varargin)
%
% Inputs
% ** startScanNr    Scan number from which to start
%
% *optional* (parameter, value) pairs
% ** ('stopScanNr', scanNum)    Final scan number, default = inf
% ** ('step', 1)                Scan number step size, default = 1
% ** ('detector','pilatus_1/2')   Detector group inside the H5 file
% ** ('outdirData','path/')     Output path default = '~/Data10/analysis/radial_integration/'
% ** ('integMask', 'path/filename') Azimuthal integration mask, default = '~/Data10/analysis/data/integration_masks.mat'
% ** ('basePath', '~/Data10')   Base path to provide to utils.find_nexus_file, default = '~/Data10'   
% ** ('integMaskPath', [])      Path to the fast integration routine, default = '<this file>/../+math/integMask/integMask'
% ** ('fastInteg', true)        Use the fast integration routine, default = true
% ** ('online', true)           Check whether the raw data files are ready. Only needed during a beamtime and should be disabled for offline
%                               processing for better performance, default = true
% ** ('useWeights', false)      Use pixel-based weights during the integration, default = false
%
% EXAMPLES:
%   beamline.integrate_range_nexus(5);
%   beamline.integrate_range_nexus(5,'detector','pilatus_2');
% 
% see also: utils.find_nexus_file, beamline.radial_integ


%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2019 by Paul Scherrer Institute (http://www.psi.ch)    |
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


function integrate_range_nexus(startScanNr, varargin)

par = inputParser;
par.addParameter('stopScanNr', inf, @isnumeric)
par.addParameter('step', 1, @isnumeric)
par.addParameter('detector', '', @ischar)
par.addParameter('outdirData', '~/Data10/analysis/radial_integration', @ischar)
par.addParameter('integMask', '~/Data10/analysis/data/integration_masks.mat', @ischar)
par.addParameter('basePath', '~/Data10', @ischar)
par.addParameter('integMaskPath', [], @ischar)
par.addParameter('fastInteg', true, @islogical)
par.addParameter('online', true, @islogical)
par.addParameter('useWeights', false, @islogical)


par.parse(varargin{:})
vars = par.Results;
vars.startScanNr = startScanNr;

if isempty(vars.detector)
    error('Specify the detector to integrate');
end
% load integration mask
[~, ~, ext] = fileparts(vars.integMask);
switch ext
    case '.h5'
        s = io.HDF.hdf5_load(vars.integMask, '-c');
        s.integ_masks = beamline.convert_integ_mask(s.integ_masks);
    case '.mat'
        s = load(vars.integMask);
    otherwise
        error('Unknown file extension.')
end

% not sure why one would need that at this point but let's keep it for now
fprintf('center at (x, y) = (%.1f, %.1f)\n',s.center_xy(1), s.center_xy(2));


scan = startScanNr;

if vars.fastInteg
    if isempty(vars.integMaskPath)
        vars.integMaskPath = fullfile(fileparts(mfilename('fullpath')), '../+math/integMask/integMask');
    end
    % try to load the first dataset to get the detectorPath
    fname = utils.find_nexus_file(vars.basePath, scan);
    while ~utils.nexus_file_is_ready(fname)
        fprintf('Waiting for scan %05d.\n', scan);
        pause(5);
        fname = utils.find_nexus_file(vars.basePath, scan);
    end
    h = io.nexus_read(fname, 'filter', vars.detector);
    
    detectorPath = fullfile(h.NXlocation, 'data');
    
    if vars.useWeights
        weights = ' -w';
    else
        weights = '';
    end

end


while scan <= vars.stopScanNr
    integTime = tic();
    fprintf('Scan S%05d: ', scan);
    if vars.online
        fprintf('Checking nexus file. ');
        fname = utils.find_nexus_file(vars.basePath, scan);
        while ~utils.nexus_file_is_ready(fname)
            fprintf('Waiting for scan %05d.\n', scan);
            pause(5);
            fname = utils.find_nexus_file(vars.basePath, scan);
        end

    else
        subdir = fullfile(vars.basePath, utils.compile_x12sa_dirname(scan));
        fname = dir(fullfile(subdir, '*.h5'));
        fname = fullfile(subdir, fname(1).name);
    end
    prepTime = toc(integTime);
    fprintf('Performing azimuthal integration. ');
    try
        if vars.fastInteg
            [~, fname_prefix] = fileparts(fname);
            outputFilename = fullfile(vars.outdirData, [fname_prefix '_00000_00000_integ.h5']);
            system(sprintf('%s -m %s -d %s -p %s -o %s %s ', vars.integMaskPath, vars.integMask, fname, detectorPath, outputFilename, weights));
        else
            beamline.radial_integ_nexus(fname, 'integMask', s, 'outdirData', vars.outdirData, 'detector', vars.detector);
        end
    catch ME
        if vars.online
            ME.rethrow();
        else
            fprintf('\n\n Azimuthal integration failed for scan %05d.\n', scan);
        end
    end
    integTimeFull = toc(integTime);
    fprintf('Finished after %.2f s (%.2f/%.2f)\n', integTimeFull, prepTime, integTimeFull-prepTime);
    scan = scan + vars.step;
end





end