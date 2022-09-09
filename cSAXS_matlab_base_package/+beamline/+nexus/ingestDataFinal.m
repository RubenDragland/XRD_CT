%INGESTDATAFINAL ingest the last dataset and miscellaneous data
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


function s = ingestDataFinal(s)
import beamline.nexus.*

if s.ingestion.ingestData
    fprintf('Ingesting last raw dataset, derived data and miscellaneous.\n');
    % make sure that the path is not empty and ends with '/'
    if ~exist(s.ingestion.tmpDir, 'dir')
        mkdir(utils.abspath(s.ingestion.tmpDir))
    end
    if isempty(s.ingestion.tmpDir)
        error('s.ingestion.tmpDir cannot be empty!');
    end
    if ~strcmp(s.ingestion.tmpDir(end), '/')
        s.ingestion.tmpDir = [s.ingestion.tmpDir '/'];
    end
       
    S = io.spec_read(s.spec.specDatFile, 'ScanNr', s.scanNr);
    if s.verbose_level<2
        fprintf(' Ingestion of dataset %u', S.dataset_id);
    end
    exclude = '';
    % prepare exclude string
    for ii=1:length(s.ingestion.exclude)
        exclude = [exclude ' --exclude ' s.ingestion.exclude{ii}];
    end
    specFiles = getSpecFiles(s.spec.specDatFile, S.dataset_id);
    if ~s.debug
        for ii=1:length(specFiles)
            [~, out] = system(sprintf('%s -o %s --ingest --online --archive --eaccount %s --dsetID %u --specFile %s %s --specRead %s', s.ingestion.cSAXSIngest, s.ingestion.tmpDir, s.eaccount, S.dataset_id, specFiles{ii}, exclude, s.spec.specParser));
        end
        [~, out] = system(sprintf('%s -o %s --ingest --online --archive --eaccount %s --skipDetector %s --specRead %s', s.ingestion.cSAXSIngest, s.ingestion.tmpDir, s.eaccount, exclude, s.spec.specParser));
    else
        fprintf('%s -o %s --ingest --online --archive --eaccount %s --dsetID %u --specFile %s %s --specRead %s\n', s.ingestion.cSAXSIngest, s.ingestion.tmpDir, s.eaccount, S.dataset_id, s.spec.specDatFile, exclude, s.spec.specParser);
        fprintf('%s -o %s --ingest --online --archive --eaccount %s --skipDetector %s --specRead %s', s.ingestion.cSAXSIngest, s.ingestion.tmpDir, s.eaccount, exclude, s.spec.specParser);
    end
    
end

end


function specFiles = getSpecFiles(specDatFile, datasetID)

res = strsplit(specDatFile, '/');
path = strjoin(res(1:end-1), '/');

specPattern = fullfile(path, '/*.dat');

cmd = sprintf('grep -l ''#C meta dataset_id int 1 %d '' %s ', datasetID, specPattern);
[~,sysout] = system(cmd);

specFiles = strsplit(sysout, '\n');
specFiles = specFiles(1:end-1);


end