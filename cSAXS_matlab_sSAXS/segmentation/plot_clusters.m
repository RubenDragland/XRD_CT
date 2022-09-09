% PLOT_CLUSTERS Extracts the representative signals from each cluster
% the representative signals are those in the cluster that are farthest away
% from other clusters.
%
% [res_segment] = plot_clusters(cluster, prep_data, par)
%
%   Inputs:
%      **cluster    indexing results of clustering analysis
%      **prep_data  input data I(q)
%      **par        structure with flags, q-range, save, and plotting parameters
%   *returns*          
%      ++res_segment    structure with matrices resulting from segmentation


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

function [res_segment] = plot_clusters(cluster, prep_data, par)
% plot clusters
clusters = cluster.num_clusters;
color_m = par.color_map;

ct = 1;
res_segment = struct();
for ii = 1:length(prep_data)
    % select size of first sample
    size_data = [prep_data(ii).scan_dimension(2), prep_data(ii).scan_dimension(1)] + 1;
    res_segment(ii).size_data = size_data;
    res_segment(ii).sample_name = prep_data(ii).sample_name;
    
    % select the range of sample
    if ii == 1
        first_idx = 1;
        last_idx = size(prep_data(ii).sym_int, 1);
    else
        first_idx = 0;
        last_idx = 0;
        for jj = 2:ii
            first_idx = first_idx + size(prep_data(jj-1).sym_int,1);
            last_idx = (first_idx) + size(prep_data(jj).sym_int,1);
            if jj == ii
                first_idx = first_idx +1;
            end
            
        end
    end
    % select the index
    idx_1 = [];
    idx_1 = cluster.idx(first_idx:last_idx);
    
    % plot in different figures
    fig1 = figure(65);
    for kk = 1:clusters
        subplot(length(prep_data), clusters, ct)
        index = find(idx_1 == kk);
        im_cluster = zeros(size_data(1)*size_data(2),1);
        im_cluster(index, 1) = 1;
        im_cluster = reshape(im_cluster, size_data(1), size_data(2));
        im_cluster = utils.adjust_projection(im_cluster, par.par_sSAXS.snake_scan, par.par_sSAXS.fast_axis_x);
        f = sprintf('cluster%d', kk);
        res_segment(ii).(f) = im_cluster;
        temp1 = zeros(size_data(1),size_data(2),3);
        % define the colors
        temp1(:,:,1) = im_cluster.*color_m(kk, 1);
        temp1(:,:,2) = im_cluster.*color_m(kk, 2);
        temp1(:,:,3) = im_cluster.*color_m(kk, 3);
        imagesc(temp1);
        box on
        axis tight equal xy
        title(sprintf('%s', prep_data(ii).sample_name), 'Interpreter', 'none')
        set(gca, 'FontSize', 8)
        axis off
        ct = ct + 1;
    end
    
    % plot all in one
    fig2 = figure(75);
    temp = zeros(size_data(1),size_data(2),3);
    for kk = 1:clusters
        index = find(idx_1 == kk);
        im_cluster = zeros(size_data(1)*size_data(2),1);
        im_cluster(index, 1) = 1;
        im_cluster = reshape(im_cluster, size_data(1), size_data(2));
        im_cluster = utils.adjust_projection(im_cluster, par.par_sSAXS.snake_scan, par.par_sSAXS.fast_axis_x);
        % define the colors
        for iii = 1:size_data(1)
            for jjj = 1:size_data(2)
                if im_cluster(iii, jjj) == 1 
                    temp(iii, jjj, 1) = color_m(kk, 1);
                    temp(iii, jjj, 2) = color_m(kk, 2);
                    temp(iii, jjj, 3) = color_m(kk, 3);
                end
            end
        end
    end
    res_segment(ii).map_color = temp;
    subplot(1, length(prep_data), ii)
    imagesc(temp);
    box on
    axis tight equal xy
    title(sprintf('%s', prep_data(ii).sample_name), 'Interpreter', 'none')
    set(gca, 'FontSize', 8)
    axis off
end

name_fig = fullfile(par.save_path, sprintf('segmented_clusters_%s_%s.jpg', par.dataset_name, par.add_name));
print(fig1, name_fig,'-djpeg','-r600');

name_fig = fullfile(par.save_path, sprintf('segmented_samples_all_%s_%s.jpg', par.dataset_name, par.add_name));
print(fig2, name_fig,'-djpeg','-r600');
