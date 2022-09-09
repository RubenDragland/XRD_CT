%APPENDEIGER500K convert raw data into h5 file
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


function s = appendEiger500K(s)
import beamline.nexus.*

load_dir = fullfile(s.eiger500K.path, s.load_dir);

if ~exist(load_dir, 'dir')
    if s.verbose_level > 3
        warning('Raw data path for Eiger 500K %s not found', load_dir)
    end
    return
end
s.rawFileExists = true;

% prepare output directory
prepareTargetDir(s);


if ~exist(s.targetFile, 'file') || ~io.HDF.hdf5_dset_exists(s.targetFile, 'eiger_3', '/entry/instrument')
    if s.verbose_level<2
        fprintf(' Eiger 500K ');
    end
    % get all file names
    file_list = dir(fullfile(load_dir, [s.eaccount '*.h5']));
    if length(file_list)==1
        % handle single H5 files separately: it is faster to just copy the
        % H5 entry instead of reading -- decompressing -- writing
        if ~s.debug
            system(sprintf('h5copy -p -i %s -o %s -s /eh5/images -d /entry/instrument/eiger_3/data', fullfile(file_list(1).folder, file_list(1).name), s.targetFile));
            if s.eiger500K.delete_raw_data
                s.deleteData.file_mask{end+1} = fullfile(file_list(1).folder, file_list(1).name);
            end
        end
    else
        num_files = length(file_list);
        tmp.entry.instrument.eiger_3.data = zeros(1030,514,num_files);
        for ii=1:num_files
            tmp.entry.instrument.eiger_3.data(:,:,ii) = io.HDF.hdf5_load(fullfile(file_list(ii).folder, file_list(ii).name), '/eh5/images');
        end
        io.HDF.save2hdf5(utils.abspath(s.targetFile), tmp, 'comp', 1);
    end
end
end