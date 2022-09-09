%APPENDEIGER1P5M convert raw data into h5 file
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


function s = appendEiger1p5M(s)
import beamline.nexus.*

% import from s
converter_path = s.eiger1p5M.converter_path;
scan = s.scanNr;

load_dir = fullfile(s.eiger1p5M.path, utils.compile_x12sa_dirname(scan));

if ~exist(load_dir, 'dir') || isempty(dir(fullfile(load_dir, '*.raw')))
    if s.verbose_level > 3
        warning('Raw data for Eiger1.5M %s not found', load_dir)
    end
    return
end
s.rawFileExists = true;

% prepare output directory
prepareTargetDir(s);




% get the hdf5 files in the output directory
if exist(s.targetFile, 'file')
    list_h5 = dir(s.targetFile);
else
    list_h5 = [];
end

if isempty(list_h5) || ~io.HDF.hdf5_dset_exists(s.targetFile, 'data', '/entry/instrument/eiger_4/')

    if s.verbose_level<2
        fprintf(' Eiger1p5M ');
    end
    % file does not exist, let's create a new one
    
    list_raw = dir(fullfile(load_dir, 'run_d0_f0000000*.raw'));
    
    
    if ~s.debug
        systemcall = [converter_path ' ' fullfile(list_raw(1).folder,list_raw(1).name) ' ' s.targetFile];
        %fprintf('%s\n',systemcall);
        [stat,out] = system(systemcall);
        
        list_h5 = dir(s.targetFile);
        if isempty(list_h5)
            error(sprintf('After conversion did not find any h5 in %s\n',load_dir))
            return
        end
        if numel(list_h5)>1
            error(sprintf('After conversion I found more than one h5 in %s\n',load_dir))
            return
        end
        
        h5fileinfo = h5info(fullfile(list_h5.folder,list_h5.name), '/entry/instrument/eiger_4/data');
        
        nframes_converted = h5fileinfo.Dataspace.Size(3);
        out = splitlines(out);
        nframes_expected = str2num(out{end-2}(14:end));
        %             fprintf('Frames expected: %i, frames converted %i \n', nframes_expected, nframes_converted)
        if nframes_converted == nframes_expected
            if s.verbose_level > 2
                fprintf('Eiger1p5M: Scan %i succefully converted to H5.\n', scan);
            end
            if s.eiger1p5M.delete_raw_data
                s.deleteData.file_mask{end+1} = sprintf('%s/*.raw',load_dir);
            end
        else
            error('Scan %i WAS NOT CONVERTED to H5\n', scan)
%             delete(sprintf('%s/*.h5',load_dir))
        end
    end
end
end