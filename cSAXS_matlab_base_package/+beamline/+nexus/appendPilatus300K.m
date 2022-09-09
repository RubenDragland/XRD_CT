%APPENDPILATUS300K read cbf data and append to already existing file or
%create a new NEXUS file. 
% 
% ** s              convert2NEXUS data structure
%
% returns:
% ++ s              updated convert2NEXUS data structure
%
% see also: beamline.nexus.convert2NEXUS

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

function s = appendPilatus300K(s)
import beamline.nexus.*

load_dir = fullfile(s.pilatus300k.path, s.load_dir);

if ~exist(load_dir, 'dir')
    if s.verbose_level > 3
        warning('Raw data path for Pilatus300K %s not found', load_dir)
    end
    return
end

s.rawFileExists = true;

% prepare output directory
prepareTargetDir(s);

if ~exist(s.targetFile, 'file') || ~io.HDF.hdf5_dset_exists(s.targetFile, 'pilatus_2', '/entry/instrument')

    if s.verbose_level<2
        fprintf(' Pilatus300K ');
    end
    res = dir(fullfile(load_dir,'*.cbf'));
    if ~isempty(res)
        % prepare target file if it does not exist
        if isempty(s.targetFile)
            s.targetFile = fullfile(s.raw_data_path, s.targetDir, [res(1).name(1:end-4) '.h5']);
        end
        if ~s.debug
            % read and write detector data
            [stat, out] = system([s.pilatus300k.converter ' -s ' load_dir ' -o ' utils.abspath(s.targetFile) ' -p /entry/instrument/pilatus_2']);
            if stat
                fprintf(out);
                error('Failed to append Pilatus300K.');
            end

            if s.pilatus300k.delete_raw_data
                s.deleteData.file_mask{end+1} = fullfile(load_dir,'*.cbf');
            end
        end
    end
end
end