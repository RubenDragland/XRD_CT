% CALCULATE_PCA estimates the principal components from the dataset
%
%   [PCA] = calculate_PCA(input, data, par)
%   Inputs:
%       **input   structure with first and second derivatives 
%       **data    input data I(q)
%       **par     structure with q-axis, save, and plotting parameters
%   *returns*
%       ++ PCA    calculated PCA

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

function [PCA] = calculate_PCA(input, data, par)

%% prepare PCA input
temp = [input.idx1];
pca_input = reshape(temp, size(par.xq, 2), size(data, 1))';

if par.use_2nd_derivative
    temp = [input.idx2];
    temp = reshape(temp, size(par.xq, 2), size(data, 1));
    pca_input = [pca_input temp'];
end

% remove zeros
pca_input(:, ~any(pca_input,1)) = [];
%% step 2.7: calculate the PCA
[~,PCA.score,~,~,PCA.explained,~] = pca(pca_input, 'VariableWeights', 'variance',  'Centered', false);%

% plot the results of PCA: three first components
figure
plot3(PCA.score(:,1), PCA.score(:,2), PCA.score(:,3), '.');
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
grid on
box on
axis tight
set(gca, 'FontSize', 20, 'LineWidth', 2)

% create the folder
if exist(par.save_path)~=7
    mkdir(par.save_path)
end

name_fig = fullfile(par.save_path, sprintf('PCA_result_%s_%s.jpg', par.dataset_name, par.add_name));
print(gcf, name_fig,'-djpeg','-r600');
%% plot the L-curve
figure
semilogy(PCA.explained, 'ob--', 'MarkerFaceColor', 'b','MarkerSize', 15)
ylabel('Proportion of variance');
xlabel('PC_i');
set(gca, 'FontSize', 25, 'LineWidth', 3)
grid on
axis tight

name_fig = fullfile(par.save_path, sprintf('PCA_variance_%s_%s.jpg', par.dataset_name, par.add_name));
print(gcf, name_fig,'-djpeg','-r600');