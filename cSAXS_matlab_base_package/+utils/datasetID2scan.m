% DATASETID2SCAN searches through a spec file and finds scan numbers
% associated with a dataset_id
% 
% ** specPath       path to the spec file or spec directory
%
% *optional*
% ** scanNr         scan number or scan number range to limit the search (-> faster)
% ** dataset_id     limit search to a specific dataset_id or range
%
% returns:
% ++ D              cell array of dataset_ids and scan numbers
%
% EXAMPLES:
%   D = utils.datasetid2scan('~/Data10');
%   D = utils.datasetid2scan('~/Data10', 'scanNr', [1:200]);
%   D = utils.datasetid2scan('~/Data10', 'dataset_id', [1:200]);
%   D = utils.datasetid2scan('~/Data10', 'scanNr', [1:200]', dataset_id', 16);
%
% see also: io.spec_read

% Academic License Agreement
%
% Source Code
%
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS 
%   ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
%
% Terms and Conditions of the LICENSE
% 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
%       hereinafter set out and until termination of this license as set forth below.
% 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
%       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
%       LICENSEE’s responsibility to ensure its proper use and the correctness of the results.”
% 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
%       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
%       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
%       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
% 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
%       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
%       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
%       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
%       in the commercial use, application or exploitation of works similar to the PROGRAM.
% 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
%       another computing language:
%       "Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
%
% Additionally, any publication using the package, or any translation of the code into another computing language should cite for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089).
% 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
%       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
% 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
%       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
%       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
%       © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
% 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
% 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
%       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
%       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
%       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
% 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
%       the courts of Zürich, Switzerland. 

function [S] = datasetID2scan(specPath, varargin)

    par = inputParser;
    par.addParameter('scanNr', -1, @isnumeric)  % scan numbers
    par.addParameter('dataset_id', -1, @isnumeric)  % dataset_ids
    
    par.parse(varargin{:})
    vars = par.Results;
        
    specDatFile = beamline.find_specDatFile(specPath);
    
    if any(vars.scanNr < 0)
        % find first scan number
        [~, out] = system(sprintf('grep ''#S'' %s | head -1', specDatFile));
        out = strsplit(out, '\n');
        out = strsplit(out{1}, ' ');
        start = str2double(out{2});
        
        startSpec = max(start, 1);
        
        % find last scan number
        [~, out] = system(sprintf('grep ''#S'' %s | tail -1', specDatFile));
        out = strsplit(out, '\n');
        out = strsplit(out{1}, ' ');
        endSpec = str2double(out{2});
        fprintf('Checking scans %u to %u.\n', startSpec, endSpec);
        vars.scanNr = [startSpec:endSpec];
    end
    
    specData = io.spec_read(specDatFile, 'ScanNr', vars.scanNr);
    S = [];
    datasetID = [];
    for ii=1:length(specData)
        if isfield(specData{ii}, 'dataset_id')
            if any(vars.dataset_id < 0) || any(ismember(vars.dataset_id, specData{ii}.dataset_id))
                if isempty(datasetID) || specData{ii}.dataset_id>datasetID
                    datasetID = specData{ii}.dataset_id;
                    S{end+1}.dataset_id = datasetID;
                    S{end}.scanNr = getScanNr(specData{ii}.S);
                else
                    S{end}.scanNr = [S{end}.scanNr getScanNr(specData{ii}.S)];
                end
            end
        end
    end
    
    
    
end


function scanNr = getScanNr(S)
scanNr = strsplit(S, ' ');
scanNr = str2double(scanNr{2});
end




















