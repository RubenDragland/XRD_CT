%STARTCONVERSION while loop; converting and appending data into nexus hdf5
%files
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

function startConversion(s)
import beamline.nexus.*

% data ingestion needs an X server
if s.ingestion.ingestData
    [~, out] = system('echo $DISPLAY');
    if isempty(strtrim(out))
        error('The data ingestion needs an X server. Please restart Matlab with display enabled.')
    end
    pgroup = strrep(s.eaccount, 'e', 'p');
    res = dir(fullfile(s.ingestion.tmpDir, pgroup, '.cscs'));
    if isempty(res)
        if ~exist('~/Documents/','dir')
            mkdir('~/Documents/')
        end
        error(sprintf('A new token for archiving with CSCS needs to be requested. Please open a new terminal and run \n\t %s -e %s -o %s --token', s.ingestion.cSAXSIngest, s.eaccount, s.ingestion.tmpDir));
    end
end

% enable ptycho-exit routine
s.status = true;

% check that first scan number exists
if s.scanNr == 1
    while true
        try
            dummy = io.spec_read(s.raw_data_path, 'ScanNr', s.scanNr, 'CReader', false);
            break;
        catch
            s.scanNr = s.scanNr + 1;
        end
    end
end

% final ingestion
if s.finished
    m=input('====== Final ingestion ======\nPlease make sure that the data processing has finished before you continue.\nDo you want to continue, Y/N [Y]:','s');
    if lower(m)=='n'
        error('Stopped by user.')
    end
end

% let's make sure that the output directory exists
if ~exist(fullfile(s.raw_data_path, s.targetDir), 'dir')
    mkdir(fullfile(s.raw_data_path, s.targetDir))
end

% ingestion only works with dataset_ids
if ~s.skipNonDatasets && s.ingestion.ingestData
    error('Cannot ingest datasets without dataset_id. Either enable skipNonDatasets or disable ingestData.')
end

% prepare default flags and storages
s.rawFileExists = false;
s.firstScanNumber = s.scanNr;
waitingStatus = false;
finishup = utils.onCleanup(@(x) converter_exit(x), s);

while true
%     tic;
    % make sure that the file masks are empty
    s.deleteData.file_mask = [];
    s.deletedData = false;

    % check if the next scan has started and find the next scan number
    [started, newScan, s.spec.specDatFile] = beamline.next_scan_started(s.raw_data_path, s.scanNr);
    
    if started
        skip = false;
        if s.skipNonDatasets
            % check dataset_id
            try
                S = io.spec_read(s.spec.specDatFile, 'ScanNr', s.scanNr);
                if ~isfield(S, 'dataset_id')
                    skip = true;
                end
                s.tmp.S = S;
            catch
                skip = true;
            end
        end
        % run function sequence to create nexus file
        if waitingStatus
            fprintf('\n');
        end
        if ~skip
            if waitingStatus
                fprintf('\n');
            end
            waitingStatus = false;
            s.load_dir = utils.compile_x12sa_dirname(s.scanNr);
            s.targetPath = fullfile(s.raw_data_path, s.targetDir, utils.compile_x12sa_dirname(s.scanNr, s.targetDirIsFlat));
            s.targetFile = fullfile(s.targetPath, sprintf(s.fnamePattern, s.scanNr));
            
            s.status = 1;
            finishup.update(s);
            fprintf('Converting scan %d: ', s.scanNr);
            for funcSeqIter=1:length(s.funcSequence)
                finishup.update(s);
                s = s.funcSequence{funcSeqIter}(s);
            end
            s.status = 0;
            finishup.update(s);
        else
            warning('No dataset_id found for scan S%05d. Skipping conversion.', s.scanNr);
        end
        s.scanNr = newScan;
        if s.verbose_level<2
            fprintf('\n');
        end
    elseif s.finished
        % wrap-up the beamtime
        if waitingStatus
            fprintf('\n');
        end
        s.load_dir = utils.compile_x12sa_dirname(s.scanNr);
        s.targetPath = fullfile(s.raw_data_path, s.targetDir, utils.compile_x12sa_dirname(s.scanNr, s.targetDirIsFlat));
        s.targetFile = fullfile(s.targetPath, sprintf(s.fnamePattern, s.scanNr));
        
        fprintf('Converting scan %d: ', s.scanNr);
        for funcSeqIter=1:length(s.funcSequence)
            finishup.update(s);
            s = s.funcSequence{funcSeqIter}(s);
        end
        finishup.update(s);
        fprintf('\n');
        s = beamline.nexus.ingestDataFinal(s);
        finishup.update(s);
        break;
        
        
    else
        % let's wait until the next scan is ready
        if ~waitingStatus
            fprintf('\nWaiting for next scan to start... ');
            waitingStatus = true;
        end
        
        fprintf('\b');
        printWaiting();
        pause(1);
    end
    s.tmp = [];
%     toc;
end

fprintf('\n\n');
fprintf('Done.\n')

s.status = 0;       % Successfully finished no need to cleanup
finishup.update(s);



end

function printWaiting()
persistent counter
if isempty(counter)
    counter = 0;
end

vals = {'|', '/', '-', '\\'};

fprintf(vals{counter+1});

counter = mod(counter+1,4);

end