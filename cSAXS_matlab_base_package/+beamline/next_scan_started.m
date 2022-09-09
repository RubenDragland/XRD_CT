%NEXT_SCAN_STARTED if the next scan has started, the first return value is
% true.
% [started, scanNr] = next_scan_started(specDatFile,scanno);
%
% ** specDatFile            path to the SPEC file / spec directory
% ** scanno                 current scan number
%
% returns:
% ++ started                true if scanno is not the last scan
% ++ nextScanno             next scan number
% ++ specDatFile            specDatFile used to determine if the next scan has started
%
% see also: beamline.find_specDatFile


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

function [started, nextScanno, specDatFileOut] = next_scan_started(dataPath,scanno)

[~, multipleSpecDatFiles] = beamline.find_specDatFile(dataPath);
nextScanno = -1;
if length(multipleSpecDatFiles)==1
    singleFile = true;
else
    singleFile = false;
end
specPattern = fullfile(multipleSpecDatFiles(1).folder, 'specES1*.dat');

cmd = sprintf('grep ''#S '' %s | grep -A 1 ''#S %d ''', specPattern, scanno);
[~,sysout] = system(cmd);
sub = strsplit(sysout, '\n');

if length(sub)<3
    arroutOld = regexp(sub{1},'[:\n ]','split');
    % check if this scan is the last scan in all files
    lastScan = getLastScan(specPattern, singleFile);
    if lastScan==scanno
        started = false;
    else
        % apparently there are more scan numbers
%         error('Scan number %d does not exist.', scanno);
        started = false;
    end
    
else
    arroutOld = regexp(sub{end-2},'[:\n ]','split');
    arrout = regexp(sub{end-1},'[:\n ]','split');
    if length(arrout)<4
        % no scan info available; check next scan
        lastScan = getLastScan(specPattern, singleFile);
        nextScanno = findNextTrueScan(specPattern, scanno, singleFile);
        
    else
        if length(sub)>5
            nextScanno = findClosestScan(sub, scanno, singleFile);
        else
            nextScanno = str2double(arrout{3-singleFile});
        end
    end
    
    if nextScanno > scanno
        started = true;
    else
        started = false;
    end
end

if singleFile
    specDatFileOut = fullfile(multipleSpecDatFiles(1).folder, multipleSpecDatFiles(1).name);
else
    specDatFileOut = arroutOld{1};
end


end

function lastScan = getLastScan(path, singleFile)

cmd = sprintf('grep ''#S '' %s | tail -1', path);
[~,sysout] = system(cmd);
sub = strsplit(sysout, '\n');
arrout = regexp(sub{1},'[:\n ]','split');
lastScan = str2double(arrout{3-singleFile});

end

function scanno = findNextTrueScan(specPattern, scanno, singleFile)

range = 2;
while true
    cmd = sprintf('grep ''#S '' %s | grep -A %d ''#S %d ''', specPattern, range, scanno);
    [~,sysout] = system(cmd);
    sub = strsplit(sysout, '\n');
    arrout = regexp(sub{end-1},'[:\n ]','split');
    if length(arrout) > 3
        scanno = str2double(arrout{3-singleFile});
        break;
    end
    range = range + 1;
end

end

function nextScan = findClosestScan(sub, scanno, singleFile)
% keyboard

scanNrs = [];
for ii=1:floor(length(sub)/3)
    arrout = regexp(sub{end-1-3*(ii-1)},'[:\n ]','split');
    scanNrs = [scanNrs str2double(arrout{3-singleFile})];
end
    
nextScan = min(scanNrs(scanNrs>scanno));

end


