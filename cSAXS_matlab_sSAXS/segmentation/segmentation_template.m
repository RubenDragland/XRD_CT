% This template is used for segmentation of scanning SAXS/WAXS images based on signal
% similarity

%% STEP 1: set the sample name and the paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.dataset_name = 'sample_name';           % same name used in SAXS_caller_template
par.base_path = '~/Data10/';                % for online: base_path = '~/Data10/';
                                            % for offline: '/das/work/p17/p17041'
par.add_name = 'ID';                        % add a name for labeling the results and images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% where the fourier components are saved
par.data_path = fullfile(par.base_path, 'analysis/fourier_components', par.dataset_name);
% where to save the results and figures
par.save_path = fullfile(par.base_path, 'analysis/segmentation', par.dataset_name);
% check the samples available
par.list = dir(par.data_path);

% show list
fprintf('***************************************************************** \n')
for ii = 1:length(par.list) % start in 3 because of dir
    fprintf('scan %d, name = %s \n', ii-2, par.list(ii).name)
end
fprintf('***************************************************************** \n')

%% STEP 2: prepare input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.select_list = [7:12];       % = [5, 6] to select specific samples
                                % = [] to select all samples in par.list
par.interpolate_log = 0;        % = 1 to interpolate logarithmically (used for SAXS)
par.interpolate_lin = 1;        % = 1 to interpolate linearly (usually used for WAXS)
par.interpolate_ratio = 10;     % = 10, reduces the points in pixel range bt 10 times
par.normalize_mean = 0;         % = 1, normalizes each I(q) by its mean value: less impact of intensity
par.pixel_range = 15:1100;      % range of pixels usually 15:1100 (exclude noisy values)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[input, prep_data, data, par] = prepare_input(par);

%% STEP 3: check the position of the index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_figures = 1;       % = 1, to show the ivatives towards 0
stepping = 5000;         % how many patterns to skip
Xaxis = 'log';          % show 'log' scale or 'lin' scale
Yaxis = 'log';          % show 'log' scale or 'lin' scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (plot_figures)
    for ii = 1:stepping:length(input)
        figure(281)
        clf
        plot(par.xq, input(ii).interpol, 'LineWidth', 0.5)
        hold on
        title(sprintf('point %d/%d', ii, length(input)))
        plot(par.xq, input(ii).idx1, '*k')
        plot(par.xq, input(ii).idx2, 'or')
        set(gca, 'YScale', Yaxis, 'XScale', Xaxis)
        grid on
        ylabel('Intensity [a.u.]')
        xlabel('Pixel range')
        legend('Data', 'First derivative', 'Second derivative')
        axis tight
        pause(0.1)
    end
end

%% STEP 4:  calculate PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.use_2nd_derivative = 1;     % =  also include the second derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[PCA] = calculate_PCA(input, data, par);

%% STEP 5: clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.cut_PCA = [];           % value found from evaluation, = [] for standard value of 0.3
par.range_clusters =  5;    % Ex: = 2 number of clusters is 2, = [2:6] to find the optimal number of
                            % clusters between 2 and 6 (takes about 30 sec to test 7 clusters/sample)
par.replicate_number = 20;  % if = 20, repeat k-mean 20 times
par.max_iteration = 100;    % if = 100, make 100 iterations for each k-mean replicate
par.color_map = [];         % if = [], uses jet values
par.plot_silhouette = 0;    % = 1 to plot silhoutte results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[cluster, par] = data_clustering(PCA, par);

%% STEP 6: representative signals and plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.reduction_factor = 10;  % = number of points to be selected: total_number_points/reduction_factor
par.furthest_pts = 1;       % = 1 to select the points furthest from the clusters centroids
par.centroid_pts = 1;       % = 1 to select the points nearest to the cluster centroids
par.y_lim = [];             % = [] auto, = [y_min, y_max] to define the limits in y    
par.offset_value = 0;       % = 0, does not add an offset to the curves, = 1 adds offset
par.Xaxis = 'log';          % show 'log' scale or 'lin' scale
par.Yaxis = 'log';          % show 'log' scale or 'lin' scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Rep_signal] = extract_rep_signals(PCA, cluster, data, par);
% plotting the results
[res_segment] = plot_clusters(cluster, prep_data, par);

%% STEP 7: save the representative signals

for ii = 1:cluster.num_clusters
    data_ii = [];
    data_ii = data((cluster.idx == ii) ,:)';
    name = fullfile(par.save_path, sprintf('points_in_cluster_%d_%s_%s.dat', ii, par.dataset_name, par.add_name));
    save(name,'data_ii', '-ASCII')
end
