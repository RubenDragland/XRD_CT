% GET_FRAMEID  Get frame IDs from SPEC metadata, in particular from
%     Dataset_ID, samplename, rotx, and/or roty.
%
% [frames] = get_frameID(base_path,varargin)
% 
% Inputs: 
%  ** base_path    e.g.  '~/Data10'
%  ** Name-value pair arguments specifying samplename, Dataset_ID, rotx and/or roty
% 
% *returns*
%  ++ frames    Frame ID numbers 
%
% EXAMPLE:
%
% % Return frames that have Dataset_ID = 45 and where the samplename is 'bone'
% frames = get_frameID(base_dir_beamtime,'Dataset_ID',45, 'samplename','bone');
%
% % Return all frames that have roty = 30
% frames = get_frameID(base_path,'roty',30);


%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
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
%   using the “cSAXS matlab package” developed by the CXS group,
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



function [frames] = get_frameID(base_path,varargin)
    if ((nargin<2)||(mod(nargin,2)~=1))
        error("Enter the correct <Entry>, <Value> pairs")
    end
    scann=[];
    for i=1:((nargin-1)/2)
        switch upper(varargin{i*2-1})
            case {'DATASET','DATASETID','DATASET_ID'}
                tofind = sprintf('#C meta dataset_id int 1 %d',varargin{i*2});
            case {'SAMPLE','SAMPLENAME','SAMPLE_NAME'}
                tofind = sprintf('#C meta samplename string 1 %s',varargin{i*2});
            case {'ROTX','ROT_X'}
                tofind = sprintf('#C meta rot_x double 1 %g',varargin{i*2});
            case {'ROTY','ROT_Y'}
                tofind = sprintf('#C meta rot_y double 1 %g',varargin{i*2});
            otherwise
        end
        if (i==1)
            scann = utils.spec_find_scans(base_path,tofind);
        else
            scann = intersect(scann,utils.spec_find_scans(base_path,tofind));
        end
    end
    if (isempty(scann))
        error("No scan found!\n");
    end
    
    tscan = io.spec_read(base_path,'ScanNr',scann);
    framelist=zeros(numel(scann),1);
    for i=1:numel(tscan)
        if isfield(tscan{i},'frame_id')
            framelist(i)=tscan{i}.frame_id;
        end
    end
    frames=unique(framelist(framelist>0));
end
