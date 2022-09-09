% EXTRACT_REP_SIGNALS extracts the purest signals, i.e. the signals within one 
% cluster that is the farthest from the other clusters
%
%   [Rep_signal] = extract_rep_signals(PCA, cluster, data, par)
%   Inputs:
%       **PCA       contains the results of PCA analysis by calculate_PCA
%       **cluster   indexing results of clustering analysis by data_clustering
%       **data      input data I(q)
%       **par       structure with flags, q-range, save, and plotting parameters
%   *returns*
%       ++Rep_signal    List of representative signals from each cluster

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

function [Rep_signal] = extract_rep_signals(PCA, cluster, data, par)

sel = PCA.score(:,1:cluster.num_comp);
clusters = cluster.num_clusters;
color_m = par.color_map;
%% select furthest points from the centroid
sumDistances = zeros(size(PCA.score, 1), cluster.num_clusters);
for ii = 1:size(PCA.score, 1)
    for iii = 1:clusters
        [temp] = sum((sel(ii,:) - cluster.centroid(iii,:)).^2, 2);
        sumDistances(ii, iii) = sumDistances(ii, iii) + temp;
    end
end
sumDistances = sum(sumDistances, 2);
%%
if (par.furthest_pts)
    % get representative signal
    M = [];
    figure
    
    Rep_signal.furthest = zeros(length(par.q), clusters);
    Rep_signal.q = par.q';
    
    off = 1;
    for ii = 1:clusters
        currentCluster = find(cluster.idx == ii);
        DistCluster = sumDistances(currentCluster);
        [DistCluster, idxDist] = sort(DistCluster, 'descend');
        range = currentCluster(idxDist(1:round(length(currentCluster)/par.reduction_factor )));
        Rep_signal.furthest(:, ii) = median(data(range, :))';
        plot(par.q, Rep_signal.furthest(:, ii)*off, 'LineWidth', 3, 'Color', color_m(ii, :))
        set(gca, 'YScale', par.Yaxis, 'XScale', par.Xaxis)
        hold on
        %axis tight
        grid on
        M = [M, {sprintf('S_%d', ii)}];
        if ii > 1
            off = off + par.offset_value*(exp(ii));
        end
    end
    legend(M, 'Location', 'north', 'Orientation','horizontal')
    set(gca, 'FontSize', 16, 'LineWidth', 1.5)
    if ~isempty(par.y_lim)
        ylim(par.y_lim)
    end
    xlabel('q [nm^{-1}]')
    ylabel('I [a.u.]')
    title('Points furthest from the centroid')
    axis tight
    name_fig = fullfile(par.save_path, sprintf('rep_signal_furthest_%s_%s.jpg', par.dataset_name, par.add_name));
    print(gcf, name_fig,'-djpeg','-r600');
end


%%
if (par.centroid_pts)
    % get representative signal
    M = [];
    figure
    
    Rep_signal.centroid = zeros(length(par.q), clusters);
    Rep_signal.q = par.q';
    
    off = 1;
    for ii = 1:clusters
        currentCluster = find(cluster.idx == ii);
        DistCluster = sumDistances(currentCluster);
        [DistCluster, idxDist] = sort(DistCluster, 'ascend');
        range = currentCluster(idxDist(1:round(length(currentCluster)/par.reduction_factor )));
        Rep_signal.centroid(:, ii) = median(data(range, :))';
        plot(par.q, Rep_signal.centroid(:, ii)*off, 'LineWidth', 3, 'Color', color_m(ii, :))
        set(gca, 'YScale', par.Yaxis, 'XScale', par.Xaxis)
        hold on
        %axis tight
        grid on
        M = [M, {sprintf('S_%d', ii)}];
        if ii > 1
            off = off + par.offset_value*(exp(ii));
        end
    end
    legend(M, 'Location', 'north', 'Orientation','horizontal')
    set(gca, 'FontSize', 16, 'LineWidth', 1.5)
    if ~isempty(par.y_lim)
        ylim(par.y_lim)
    end
    xlabel('q [nm^{-1}]')
    ylabel('I [a.u.]')
    title('Points nearest to the centroid')
    axis tight
    name_fig = fullfile(par.save_path, sprintf('rep_signal_centroid_%s_%s.jpg', par.dataset_name, par.add_name));
    print(gcf, name_fig,'-djpeg','-r600');
end

