% PREPARE_INPUT prepares the data for PCA analysis
%
% [input, prep_data, data, par] = prepare_input(par)
%
%   Inputs:
%      **par    stucture defining the data paths, flags, q-range, save, and 
%               plotting parameters
%   *returns*
%      ++input      structure with first and second derivatives to be used as input 
%                   for PCA analysis
%      ++prep_data  I(q) evaluated at the inflection points
%      ++data       original I(q) data
%      ++par        updated parameter structure


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
%
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing la this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS scanning SAXS package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% Additionally, any publication using the package, or any translation of the 
%     code into another computing language should cite:
%    O. Bunk, M. Bech, T. H. Jensen, R. Feidenhans'l, T. Binderup, A. Menzel 
%    and F Pfeiffer, “Multimodal x-ray scatter imaging,” New J. Phys. 11,
%    123016 (2009). (doi: 10.1088/1367-2630/11/12/123016)
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [input, prep_data, data, par] = prepare_input(par)

list = par.list;

if isfield(par, 'select_list') && ~isempty(par.select_list)
    range = par.select_list+2; % the 2 is because of dir
else
    range = 3:length(list);
end


% prepare data
prep_data = struct();
ct = 1;
for ii =  range
    % compile the data directory from data-path and file-name and look for
    % fourier.mat files
    temp_dir_mask = fullfile(par.data_path, list(ii).name,'*fourier.mat'); 
    temp = dir(temp_dir_mask);
    % check if none or more than one has been found 
    if (isempty(temp))
        exit('No files found matching %s',temp_dir_mask);
    end
    if (numel(temp) > 1)
        fprintf('Warning: more than one fourier.mat file found matching %s\n');
        fprintf('Chosing the first one:\n');
        for temp_ind=1:numel(temp)
            fprintf('%2d. %s\n',temp_ind, temp(temp_ind).name);
        end
        % wait briefly before this message scrolls away
        pause(3);
    end
    % load the first file (usually only one should exist, as re-processing
    % always overwrites the previous one, but earlier versions of the
    % segmentation-preparation script used different file-naming and thus
    % duplicates may occur in historic data)
    filename = fullfile(temp(1).folder, temp(1).name);
    fprintf('loading %s\n',filename);
    temp = load(filename);
    data = temp.data_process.sym_data;
    data = reshape(data,size(data,1)*size(data,2), []);
    data(:,~any(data,1)) = [];
    prep_data(ct).sym_int = data(:, par.pixel_range);
    prep_data(ct).scan_dimension = [temp.four_dat.par.ind_max_x, temp.four_dat.par.ind_max_y];
    prep_data(ct).sample_name = list(ii).name;
    ct = ct+1;
end

% save the q-range
par.q = temp.four_dat.integ.q(par.pixel_range);
par.par_sSAXS = temp.four_dat.par;

% prepare the data input
data = [];
for ii = 1:length(prep_data)
    data = [data; prep_data(ii).sym_int];
end



%% normalize and reduce the input data
input = struct();
% interpolation vector: reduce the number of pixels by 5
if (par.interpolate_log)
    par.xq = logspace(1, log10(size(data, 2)), round(size(data, 2)/par.interpolate_ratio));
elseif (par.interpolate_lin)
    par.xq = linspace(1, (size(data, 2)), round(size(data, 2)/par.interpolate_ratio));
elseif (par.interpolate_log) && (par.interpolate_lin)
    par.xq = [1:size(data, 2)];
end
    
parfor ii = 1:size(data, 1)
    % interpolate
    temp = smooth(data(ii, :)); % smooth(Y) uses the moving average method with span 5 and X=1:length(Y)
    temp = interp1(1:length(temp), temp, par.xq);
    if par.normalize_mean
        % normalize
        temp = temp/mean(temp);
        input(ii).interpol = temp;
    else
        input(ii).interpol = temp;
    end
    input(ii).first = diff(temp);
    input(ii).second = diff(diff(temp));
    fprintf('Derivating: point %d/%d \n', ii, length(data));
end

%% extract the features
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
for ii = 1:size(data, 1)
    temp = zci(input(ii).first);
    temp = temp + 1;
    input(ii).idx1 = zeros(1, size(input(ii).interpol, 2));
    input(ii).idx1(temp) = input(ii).interpol(temp);
    temp = zci(input(ii).second);
    temp = temp + 2;
    input(ii).idx2 = zeros(1, size(input(ii).interpol, 2));
    input(ii).idx2(temp) = input(ii).interpol(temp);
    fprintf('Extracting zeros: point %d/%d \n', ii, length(data));
end
