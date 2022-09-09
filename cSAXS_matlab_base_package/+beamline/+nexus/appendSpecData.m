%APPENDSPEC read and append spec, mcs and orchestra data to a nexus file
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

function s = appendSpecData(s)
import beamline.nexus.*

if  s.rawFileExists || s.spec.appendEmpty
    
    % prepare output directory
    prepareTargetDir(s);
    
    
    % prepare json stream
    dat = {};
    if exist(s.orchestra.path, 'dir')
        j = [];
        j.path = utils.abspath(fullfile(s.orchestra.path, sprintf('scan_%05d.dat', s.scanNr)));
        j.type = 'orchestra';
        j.groupName = 'orchestra';
        dat{length(dat)+1} = j;
    end
    if exist(s.mcs.path, 'dir')
        j = [];
        j.path = utils.abspath(fullfile(s.mcs.path, utils.compile_x12sa_dirname(s.scanNr), sprintf('%s_%05d.dat', s.eaccount, s.scanNr)));
        j.type = 'mcs';
        j.entry = 'mcs_data';
        j.groupName = 'mcs';
        dat{length(dat)+1} = j;
    end
    if exist(s.sgalil.path, 'dir')
        j = [];
        j.path = utils.abspath(fullfile(s.sgalil.path, sprintf('S%05d.dat', s.scanNr)));
        j.type = 'orchestra';
        j.groupName = 'sgalil';
        dat{length(dat)+1} = j;
    end
    
    if ~isempty(dat)
        if length(dat)==1
            datFiles = ['--datFiles ''' jsonencode(dat{1}) ''''];
        else
            datFiles = ['--datFiles ''' jsonencode(dat) ''''];
        end
    else
        datFiles ='';
    end
    
    
    % let's go!
    if ~exist(s.targetFile, 'file') || ~io.HDF.hdf5_dset_exists(s.targetFile, 'specES1', '/entry/collection')
        if s.verbose_level<2
            fprintf(' SPEC ');
        end
        systemcall = sprintf('%s -s %s -v 5 --scanNr %u --hdf5 --xmlLayout %s -o %s %s', utils.abspath(s.spec.specParser), utils.abspath(s.spec.specDatFile), s.scanNr, utils.abspath(s.spec.xmlLayoutFile), utils.abspath(s.targetFile), datFiles);
        if ~s.debug
            [stat, out] = system(systemcall);
            if stat~=0
                pause(5);
                [stat, out] = system(systemcall);
            end
        else
            fprintf([systemcall '\n']);
        end
    end
end
end