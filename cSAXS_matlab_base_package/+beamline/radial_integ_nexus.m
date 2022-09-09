%RADIAL_INTEG_NEXUS perform azimuthal integration on a nexus raw data file
%
% ** fname                  full path of the nexus raw data file
%
% *optional*
% ** detector               detector name (default: pilatus_1)
% ** outdirData             save path for the integration file
% ** integMask              full path of the integration mask
% ** rMaxForced             limit the radial integration range       
% ** adjustOrientation      apply the detector orientation before
%                           performing the azimuthal integration
%
% EXAMPLE:
%   beamline.radial_integ_nexus('~/Data10/data/S00000-00999/S00001/data.h5', 'detector', 'pilatus_1') 
%
% see also: beamline.integrate_range_nexus.m


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



function radial_integ_nexus(fname, varargin)

par = inputParser;
par.addParameter('detector', 'pilatus_1', @ischar)
par.addParameter('outdirData', '', @ischar)
par.addParameter('integMask', [], @isstruct)
par.addParameter('basePath', '~/Data10', @ischar)
par.addParameter('rMaxForced', 0, @isnumeric)
par.addParameter('saveFormat', 'h5', @ischar)
par.addParameter('adjustOrientation', false, @islogical)

par.parse(varargin{:});
vars = par.Results;

% make sure the saveFormat does not start with '.'
if strcmpi(vars.saveFormat(1), '.')
    vars.saveFormat = vars.saveFormat(2:end);
end

s = vars.integMask;

% limit radial range
if (vars.rMaxForced > 0)
    ind = find( s.integ_masks.radius < vars.rMaxForced );
    if (length(ind) < 1)
        fprintf('No radii below rMaxForced = %d found\n',vars.rMaxForced);
        return;
    end
    s.integ_masks.radius = s.integ_masks.radius(1:ind(end));
    s.integ_masks.norm_sum = s.integ_masks.norm_sum(1:ind(end), :);
end

% check if radius exists; if not, use q
if isfield(s.integ_masks,'radius')
    radius = s.integ_masks.radius;
    ind_r_max = length(radius);
else
    radius = [];
    q = s.integ_masks.q;
    ind_r_max = length(q);
end
if isfield(s.integ_masks, 'q')
    q = s.integ_masks.q;
else
    q = [];
end


% load raw data
if vars.adjustOrientation
    data = io.nexus_read(fname, 'filter', vars.detector);
else
    data = io.nexus_read(fname, 'filter', vars.detector, 'orientation', [0 0]);
end

if isempty(data)
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AZIMUTHAL INTEGRATION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[s.I_all, s.I_std] = radial_integ_mex(int32(ind_r_max),...
    int32(s.no_of_segments), s.integ_masks.norm_sum, ...
    s.integ_masks.indices, double(data.data));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isempty(vars.outdirData))
    % add slash to output directory
    if (vars.outdirData(end) ~= '/')
        vars.outdirData = [ vars.outdirData '/' ];
    end
    
    % create output directory
    [mkdir_stat,mkdir_message] = mkdir(vars.outdirData);
    if (~mkdir_stat)
        error('invalid directory %s: %s',vars.outdirData,mkdir_message);
    end
%     if ((mkdir_stat) && (isempty(mkdir_message)))
%         fprintf('The output directory %s has been created.\n',vars.outdirData);
%     else
%         fprintf('The output directory is %s.\n',vars.outdirData);
%     end
else
    fprintf('data are not saved\n');
end

% save data, if this option is enabled
if (~isempty(vars.outdirData))
    if  isfield(s, 'I_all')
        % use first file as file-name base
        [~, name] = fileparts(fname);

        fname_out = fullfile(vars.outdirData, [ name '_00000_00000_integ.' vars.saveFormat]);
%         fprintf('saving %s\n',fname_out);
        
        % remove directory information before storing the filenames
        s.integ_masks = rmfield(s.integ_masks, 'indices');
        s.integ_masks = rmfield(s.integ_masks, 'weights');
        s = rmfield(s, 'asize');
        switch vars.saveFormat
            case 'mat'
                save(fname_out,'s');
            case 'h5'
                io.HDF.save2hdf5(fname_out, s, 'overwrite', 1);
            otherwise
                error('Unknown file extension %s.', vars.saveFormat)
        end
    else
        fprintf('No data to save for directory %s\n',data_dir);
    end
end

end