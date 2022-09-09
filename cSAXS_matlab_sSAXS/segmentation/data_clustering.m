% DATA_CLUSTERING apply clustering to the PCA results
% 
%   [cluster, par] = data_clustering(PCA, par)
%   Inputs:
%       **PCA   contains the results of PCA analysis by calculate_PCA
%       **par   structure with PCA threshold, save, and plotting parameters
%   *returns*
%       ++cluster   indexing results of clustering analysis
%       ++par       updated structure

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

function [cluster, par] = data_clustering(PCA, par)

if isfield(par, 'cut_PCA') && ~isempty(par.cut_PCA)
    cut_PCA = par.cut_PCA;
else
    cut_PCA = 0.3;
end

% calculate the optimal number of clusters
cluster.num_comp = length(find(PCA.explained > cut_PCA));
sel = PCA.score(:,1:cluster.num_comp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edit here:
distance = 'sqeuclidean';
range_clusters = par.range_clusters;
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%
if length(range_clusters) > 1
    rng(1);%for reproducibility
    tic
    eva = evalclusters(sel,'kmeans','silhouette', 'Klist', range_clusters, 'Distance',distance);
    toc
    clusters = eva.OptimalK;
    figure (185)
    plot(range_clusters, eva.CriterionValues, 'ob--', 'MarkerFaceColor', 'b','MarkerSize', 15)
    ylabel('Evaluation [-1 to 1]');
    xlabel('Number of clusters');
    set(gca, 'FontSize', 25, 'LineWidth', 3)
    grid on
    name_fig = fullfile(par.save_path, sprintf('eval_%s_%s.jpg', par.dataset_name, par.add_name));
    print(gcf, name_fig,'-djpeg','-r600');
else
    clusters = range_clusters;
end
cluster.num_clusters = clusters;
%% clustering by k-means
opts = statset('Display','final');
rng(1);%for reproducibility
[cluster.idx,cluster.centroid] = kmeans(sel, clusters,'Distance',distance,...
    'Replicates',par.replicate_number,'MaxIter',par.max_iteration, 'Options',opts, 'Start', 'plus');

%% plot the silhouette
if (par.plot_silhouette)
    figure
    tic
    [silh, h] = silhouette(sel, cluster.idx, distance);
    toc
    h = gca;
    h.Children.EdgeColor = [.8 .8 1];
    xlabel 'Silhouette Value'
    ylabel 'Cluster Number'
    temp = mean(silh);
    title(sprintf('Mean Silhouette %.2f', temp));
    name_fig = fullfile(par.save_path, sprintf('cluster_silhouette_%s_%s.jpg', par.dataset_name, par.add_name));
    print(gcf, name_fig,'-djpeg','-r600');
end

%% PLOT THE CLUSTERS
figure
if isfield(par,'color_map') && isempty(par.color_map)
    c = jet;
    color_m = round(linspace(1, length(c), clusters));
    color_m = c(color_m,:);
    par.color_map = color_m;
else
    color_m = par.color_map;
end
for ii = 1:clusters
    plot3(PCA.score(cluster.idx == ii, 1), PCA.score(cluster.idx == ii, 3), PCA.score(cluster.idx == ii, 2), '.', 'Color', color_m(ii,:))
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    hold on
    grid on
    box on
    set(gca, 'FontSize', 18, 'LineWidth', 2)
    set(gcf,'color','w');
    axis tight
end
name_fig = fullfile(par.save_path, (sprintf('cluster_colors_%s_%s.jpg', par.dataset_name, par.add_name)));
print(gcf, name_fig,'-djpeg','-r600');

