% comparison of the result from optimization and measured data
% compare the 2D projections

close all;clear all
%addpath(genpath('/home/fredrimu/SASTT/cxs_software/'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%base_path = '/home/fredrimu/SASTT/';%'~/Data10';         % = ~/Data10 for online analysis, provide the path for offline analysis Ex: '/das/work/p16/p16649/'
%sample_name = 's3t2_SAXS_7sect_nocart';%'sample_name';    % name given in the saxs_caller_template
%add_name = 'HA_angles';%'ID';         % additional name the optimizations: = [ ] if not needed

parent = cd ;
base_path = '/Data sets/' ; 
base_path = [parent base_path] ; % SAFE WITH FULL PATH
sample_name = 'SASTT_carbon_knot_aligned_ASTRA_corrected' ;
add_name = '';%'ID';         % additional name the optimizations: = [ ] if not needed

checkAfterAngleOpt = 0;
use360degDetector = 0;

% 
%  sample_name_all = {'projection_a_1_or_1','projection_a_1_or_2','projection_a_1_or_3',...,
%                     'projection_a_2_or_1','projection_a_2_or_2','projection_a_2_or_3','projection_a_3_or_1','projection_a_3_or_2','projection_a_3_or_3',...,
%                     'projection_a_4_or_1','projection_a_4_or_2','projection_a_4_or_3'};
% all_reg_coeffs = [0.016 0.017 0.016 0.016 0.017 0.016 0.016 0.018 0 0.009 0.040 0.015];
% 
% for nn = 10:12%:length(sample_name_all)
    
%     if nn > 1
% %         clear a
% %         clear numOfvoxels
% %         clear numOfCoeffs
%         clear p
%         clear s
%         clear projection
%     end
    
%     sample_name = 's3t2_ha002'; %sample_name_all{nn};%'sample_name';
    
    pix_size = 0.050; %set pix size for output projection. FKM
    isShale = 0; %if shale sample.
    dualData = 0;  
    
    % add_name = 'all_it_single_axis';%'ID';         % additional name the optimizations: = [ ] if not needed
    reg_coeff = 20;% all_reg_coeffs(nn); %FKM load the optimization with this reg.coeff
    
    include_IRTT = 0;               % = 1 to include the results from IRTT
    % = 0 to not include IRTT
    theta_det = pi/2;               % Ewald sphere correction, for SAXS pi/2
    which_projections = [21]; % which projections to plot, example [1:10:300]
    plot_corrected_data = 1;        % =0, if 2D projections don't have corrections
    % =1, if corrections have been applied during alignment
    apply_mask = 0;                 % = 0 does not apply mask
    % = 1 apply mask
    std_scale = [];                 % = [] does not change  the std scale
    % = 3, changes to scale to 3*std of the data (as in SAXS_caller_template)
    
    if isShale
        cw_size = [0.05];                   % = [], color wheel will be same size (about 0.16)
    else
        cw_size = [];
    end
    % = 0.09, for decreasing  the size of  the
    % color wheel
    save_video = 0;                 % = 1, to save video, = 0 does not save
    save_figure = 1;                % = 1, to save figure, = 0 does not save
    save_I_phi_plot = 1;                % FKM saves a I_phi plot at point
    I_phi_point = [14 14];              %point  
    scaleOrientation = 0;
    customColorScale = [0 1];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    addpath(fullfile(base_path, 'cxs_software','scanning_saxs'));
    % load aligned projections
    %filename = fullfile(base_path,sprintf('analysis/SASTT/%s/projection_data/SASTT_%s_aligned_ASTRA.mat', ...
    %    sample_name, sample_name));
    filename = sprintf('%s%s.mat', base_path, sample_name); %RSD: CURRENTLY WHAT IS NEEDED. DATA SETS LOCATED OUTSIDE OF PACKAGE. 
    load(filename);
    
    
    % prepare video
    if save_video
        videoname = sprintf('video_%s_%s_reg_coeff_%0.3f.avi', sample_name, add_name,reg_coeff);
        videoname = fullfile(base_path,  sprintf('analysis/SASTT/%s/SH/%s/figures/%s', sample_name, add_name, videoname));
        if ~(exist('videoname', 'file'))
            v = VideoWriter(videoname);
            v.FrameRate = 2;
            open(v)
        end
    end
    screen_size = get(0,'screensize');
    four_dat.fmt.fig_pos_sym_int =  [150 150 screen_size(3)-300 screen_size(4)-300 ];
    four_dat.fmt.fig_pos_comb_all =  [150 150 screen_size(3)-300 screen_size(4)-300 ];
    
    

    
    
    
    %% Fix shale
    
    if isShale
        for ii = 1:length(projection)                        % FOR SHALE
            projection(ii).par.ind_max_x = 71;
            projection(ii).par.ind_max_y = 28;
        end
    end
    
    %% loop over THE projections
    for kk = 1:length(which_projections)
        % select one projection
        one_projection = projection(which_projections(kk));
        one_projection.integ.theta_det = theta_det;
        
        if use360degDetector
            one_projection.data = 0.5*(one_projection.data(:,:,1:16) + one_projection.data(:,:,17:32));
        end
        
        
        
        %comparing
        fprintf('Comparing projection %05d at rotX %.1f and rotY %.1f \n', which_projections(kk), one_projection.rot_x, one_projection.rot_y);
        fprintf('** \n');
        
        if ~(include_IRTT)
            % plot the original data
            % load the measured 2D projections
            
            % Commented out FKM
            %         filename = sprintf('%sanalysis/fourier_components/%s/%s%05d_fourier.mat', base_path, ...
            %             sample_name, one_projection.fnames.rel_dir, one_projection.fnames.first_scan_no);
            %         load(filename)
            
            % if needed, change  the saturation of the images
            if ~(isempty(std_scale))
                four_dat.fmt.std_scale = std_scale;
            end
            % if needed, change  the size of the color wheel
            if ~(isempty(cw_size))
                four_dat.fmt.cw_size = cw_size;
            end
            
            % Edit FKMUERER: Hack when four_dat is missing.
            four_dat.fnames = '';
            
            
            
            
            four_dat.par = projection(1).par;
            four_dat.par_phi_det = projection(1).integ.phi_det;
            four_dat.par.q = 100; %dummy
            four_dat.fmt.interpol = 0;
            if dualData
                four_dat.f1_amp = squeeze(mean(one_projection.data_peak1,3));
                foo = fft(one_projection.data_peak1,[],3);
            else
                if isfield(one_projection,'data_peak1') & ~isfield(one_projection,'data')
                    one_projection.data = one_projection.data_peak1;
                end
                four_dat.f1_amp = squeeze(mean(one_projection.data,3));
                foo = fft(one_projection.data,[],3);
            end
            
            four_dat.f2_amp = 0; %squeeze(abs(foo(:,:,3)));
            four_dat.f2_phase = 0;%angle(foo(:,:,2));
            four_dat.degree_orient = 0; %abs(foo(:,:,3)).^2./(abs(foo(:,:,1)).^2 + abs(foo(:,:,3)).^2);
            four_dat.orientation = 0; %phase2orientation(four_dat.f2_phase,four_dat.par_phi_det,0);
            four_dat.fmt.std_scale = [1];
            four_dat.integ.phi_det = projection(1).integ.phi_det;
            four_dat.par.x_scale = pix_size;
            four_dat.par.y_scale = pix_size;
            
            % needed to only plot the data
            processing.load_data = 0; % 1 to force re-load and re-process
            processing.save_fig  = 0;
            processing.plot_data = 1;
            processing.print_fig = 0;
            processing.movie_fig = 0;
            % define which q-range to plot
            %         processing.which_qrange_to_plot = four_dat.par.which_qrange;
            
            processing.which_range_to_plot = 100; %EDIT FKM
            
            % only plot the important plots
            four_dat.fmt.fig_sym_int = 1;        % symmetric intensity over all segments: f1_amp
            four_dat.fmt.fig_asym_sym_int = 0;   % combination of orientation (colorwheel), degree of orientation (hue)
            four_dat.fmt.fig_asym_int = 0;       % amplitude of the anisotropy: f2_amp
            four_dat.fmt.fig_orientation = 0;    % angle of orientation (in degrees): f2_phase
            four_dat.fmt.fig_orient_degree = 0;  % degree of orientation from f2_amp/f1_amp
            four_dat.fmt.fig_orient_degree_histogram = 0; % histogram of the degree of orientation vs. orientation
            four_dat.fmt.fig_scatt_asym = 0;     % scattering asymmetry
            four_dat.fmt.fig_cos_dev = 0;        % how the fit of azimuthal integration deviates from a cosine
            four_dat.fmt.fig_I_pow = 0;          % slopes of I(q) in the selected q-region
            four_dat.fmt.fig_peak_pos = 0;       % peak position and amplitude
            % define positions
            
            
            if scaleOrientation
                four_dat.fmt.scale_orientation = customColorScale;
            end
            
            
            % turn title off
            four_dat.fmt.title_on = 0;
            
            if plot_corrected_data
                processing.which_qrange_to_plot = 1;
                % calculate the fourier coefficients
                if dualData
                    [four_dat] = calculate_four_dat(four_dat, one_projection.data_peak1);
                    
                else
                    [four_dat] = calculate_four_dat(four_dat, one_projection.data);
                end
            end
            
            % if we need to apply the mask
            if (apply_mask)
                if ~(plot_corrected_data)
                    % adjust sizes, if mask and original data have different
                    % sizes
                    temp = [];
                    names = fieldnames(four_dat);
                    for pp = 2:9 % only the fields needed here
                        % y dimension
                        while size(four_dat.(names{pp}), 1) < size(one_projection.window_mask, 1)
                            temp = [four_dat.(names{pp})(1, :, :); four_dat.(names{pp})];
                            four_dat.(names{pp}) = temp;
                        end
                        % x dimension
                        while size(four_dat.(names{pp}), 2) < size(one_projection.window_mask, 2)
                            temp = [four_dat.(names{pp})(:, 1, :), four_dat.(names{pp})];
                            four_dat.(names{pp}) = temp;
                        end
                    end
                end
                % compile return structure for plotting and mask
                four_dat.f1_amp = four_dat.f1_amp.*one_projection.window_mask;
                four_dat.f2_amp = four_dat.f2_amp.*one_projection.window_mask;
                four_dat.f2_phase = four_dat.f2_phase.*one_projection.window_mask;
                four_dat.orientation = four_dat.orientation.*one_projection.window_mask;
                four_dat.degree_orient = four_dat.degree_orient.*one_projection.window_mask;
            end
            
            four_dat.fmt.parallel_plotting = false;
            
            
            % plot the mesured symmetric intensity
            figh = figure(100);  set(figh,'Position',[2026 155 1636 900]); %FKM 
            clf(figh)
            axesh = gca;
            subplot(2,2,1)
            imgh = imagesc(1);
            four_dat.fmt.figh = figh;
            four_dat.fmt.imgh = imgh;
            % plot
            plot_fourier_result(four_dat, processing);
            title(sprintf('Measured sym. int., \\beta = %.1f\\circ, \\alpha =  %.1f\\circ', one_projection.rot_x, one_projection.rot_y), 'FontSize', 15);
            caxis auto %FKM
            
            % plot the measured orientation
            four_dat.fmt.fig_sym_int = 0;        % symmetric intensity over all segments: f1_amp
            four_dat.fmt.fig_asym_sym_int = 1;   % combination of orientation (colorwheel), degree of orientation (hue)
            figh = figure(100);
            subplot(2,2,3)
            imgh = imagesc(1);
            axesh = gca;
            four_dat.fmt.figh = figh;
            four_dat.fmt.imgh = imgh;
            % plot
            plot_fourier_result(four_dat, processing);
            title(axesh, sprintf('Measured orient., \\beta = %.1f\\circ, \\alpha =  %.1f\\circ', one_projection.rot_x, one_projection.rot_y), 'FontSize', 15);
            
            % Calculate 2D projections from the optimized data with SH
            
            if checkAfterAngleOpt
                filename = fullfile(base_path, sprintf('analysis/SASTT/%s/SH/%s/optimization_output/', sample_name, add_name));
                filename = fullfile(filename, sprintf('result_%s_q%d-%d_angles_%s.mat', sample_name, projection(1).par.r_sum{1}(1), ...
                projection(1).par.r_sum{1}(end), add_name));
            else
                
            filename = fullfile(base_path, sprintf('analysis/SASTT/%s/SH/%s/results/', sample_name, add_name));
            filename = fullfile(filename, sprintf('result_%s_q%d-%d_all_again_%s_rcoeff_%0.3f.mat', sample_name, projection(1).par.r_sum{1}(1), ...
                projection(1).par.r_sum{1}(end), add_name,reg_coeff));
            
            
            end
            
            
            if ~isfile(filename)
                filename = fullfile(base_path, sprintf('analysis/SASTT/%s/SH/%s/results/', sample_name, add_name));
                %filename = fullfile(filename, sprintf('result_%s_q%d-%d_all_again_%s_reg_coeff_%0.3f.mat', sample_name, projection(1).par.r_sum{1}(1),projection(1).par.r_sum{1}(end), add_name,reg_coeff));
                filename = fullfile(filename, sprintf('result_%s_q%d-%d_all_again_%s.mat', sample_name, projection(1).par.r_sum{1}(1),projection(1).par.r_sum{1}(end), add_name)); % RSD: THE CURRENT RESULT LACKS THE LAST PART OF THE FILE NAME. 
            end
            
            load(filename)
            
            
            if dualData
                tomoSz = [74 72 72]; %s3t2
                %             tomoSz = [50 50 50]; % phantom
                numOfvoxels = prod(tomoSz);
                s.a(1).data = reshape(a_out1_peak1(1:numOfvoxels),tomoSz);
                s.a(1).l = 0;
                s.a(1).m = 0;
                s.a(2).data = reshape(a_out1_peak1((numOfvoxels+1):(2*numOfvoxels)),tomoSz);
                s.a(2).l = 2;
                s.a(2).m = 0;
                s.a(3).data = reshape(a_out1_peak1((2*numOfvoxels+1):(3*numOfvoxels)),tomoSz);
                s.a(3).l = 4;
                s.a(3).m = 0;
                s.a(4).data = reshape(a_out1_peak1((3*numOfvoxels+1):end),tomoSz);
                s.a(4).l = 6;
                s.a(4).m = 0;
                s.theta.data = reshape(theta_out1,tomoSz);
                s.phi.data = reshape(phi_out1,tomoSz);
                s.mask3D = params.mask3D;
                
                
                
                
                p = params;
                p.nx = tomoSz(2);
                p.ny = tomoSz(1);
                p.nz = tomoSz(3);
                p.opt_coeff = [0 1 1 1];
                
            end
            
            % For more efficiency the loop over projections should be done
            % inside SAXS_tomo_3D_err_metric
            
            
            if isfield(one_projection,'data_peak1')
                one_projection.data = one_projection.data_peak1;
            end
            
            if ~isfield(one_projection,'window_mask')
                if dualData
                one_projection.window_mask = ones(size(one_projection.data_peak1,1),size(one_projection.data_peak1,2));
                else
                    one_projection.window_mask = ones(size(one_projection.data,1),size(one_projection.data,2));
                end
            end
            
            
            if use360degDetector
                p.phi_det = p.phi_det(1:16);
                p.numsegments = 16;
            end
            
            [~,~,proj_SH] = optimization.SAXS_tomo_3D_err_metric([],p, s, one_projection);
            % load data as the variables must be reprocessed
            processing.load_data = 1; % 1 to force re-load and re-process
            processing.plot_data = 1;
            % define which q-range to plot
            processing.which_qrange_to_plot = 1;
            
%             if use360degDetector
%                 proj_SH.projection = 0.5*(proj_SH.projection(:,:,1:16) + proj_SH.projection(:,:,17:32));
%             end
%                 
            
            
            
            % calculate the fourier coefficients
            [four_dat] = calculate_four_dat(four_dat, proj_SH.projection);
            
            
            
            
            % if needed to apply the mask
            if (apply_mask)
                % compile return structure for plotting and mask
                four_dat.f1_amp = four_dat.f1_amp.*one_projection.window_mask;
                four_dat.f2_amp = four_dat.f2_amp.*one_projection.window_mask;
                four_dat.f2_phase = four_dat.f2_phase.*one_projection.window_mask;
                four_dat.orientation = four_dat.orientation.*one_projection.window_mask;
                four_dat.degree_orient = four_dat.degree_orient.*one_projection.window_mask;
            end
            
            % plot the simulated symmetric intensity
            four_dat.fmt.fig_sym_int = 1;        % symmetric intensity over all segments: f1_amp
            four_dat.fmt.fig_asym_sym_int = 0;   % combination of orientation (colorwheel), degree of orientation (hue)
            figh = figure(100);
            subplot(2,2,2)
            imgh = imagesc(1);
            four_dat.fmt.figh = figh;
            four_dat.fmt.imgh = imgh;
            % plot
            plot_fourier_result(four_dat, processing);
            title(sprintf('Simulated sym. int., \\beta = %.1f\\circ, \\alpha =  %.1f\\circ', one_projection.rot_x, one_projection.rot_y) , 'FontSize', 15);
            caxis auto %FKM
            
            
            % plot the simulated orientation
            four_dat.fmt.fig_sym_int = 0;        % symmetric intensity over all segments: f1_amp
            four_dat.fmt.fig_asym_sym_int = 1;   % combination of orientation (colorwheel), degree of orientation (hue)
            figh = figure(100);
            subplot(2,2,4)
            imgh = imagesc(1);
            axesh = gca;
            
            four_dat.fmt.figh = figh;
            four_dat.fmt.imgh = imgh;
            % plot
            plot_fourier_result(four_dat, processing);
            title(axesh, sprintf('Simulated orient., \\beta = %.1f\\circ, \\alpha =  %.1f\\circ', one_projection.rot_x, one_projection.rot_y), 'FontSize', 15);
        else %%NOT RELEVANT BELOW THIS, ONLY FOR IRTT
            %% plot the original data
            % load the measured 2D projections
            filename = sprintf('%sanalysis/fourier_components/%s/%s%05d_fourier.mat', base_path, ...
                sample_name, one_projection.fnames.rel_dir, one_projection.fnames.first_scan_no);
            load(filename)
            
            % if needed, change  the saturation of the images
            if ~(isempty(std_scale))
                four_dat.fmt.std_scale = std_scale;
            end
            % if needed, change  the size of the color wheel
            if ~(isempty(cw_size))
                four_dat.fmt.cw_size = cw_size;
            end
            
            % needed to only plot the data
            processing.load_data = 0; % 1 to force re-load and re-process
            processing.save_fig  = 0;
            processing.plot_data = 1;
            processing.print_fig = 0;
            processing.movie_fig = 0;
            % define which q-range to plot
            processing.which_qrange_to_plot = four_dat.par.which_qrange;
            
            % only plot the important plots
            four_dat.fmt.fig_sym_int = 1;        % symmetric intensity over all segments: f1_amp
            four_dat.fmt.fig_asym_sym_int = 0;   % combination of orientation (colorwheel), degree of orientation (hue)
            four_dat.fmt.fig_asym_int = 0;       % amplitude of the anisotropy: f2_amp
            four_dat.fmt.fig_orientation = 0;    % angle of orientation (in degrees): f2_phase
            four_dat.fmt.fig_orient_degree = 0;  % degree of orientation from f2_amp/f1_amp
            four_dat.fmt.fig_orient_degree_histogram = 0; % histogram of the degree of orientation vs. orientation
            four_dat.fmt.fig_scatt_asym = 0;     % scattering asymmetry
            four_dat.fmt.fig_cos_dev = 0;        % how the fit of azimuthal integration deviates from a cosine
            four_dat.fmt.fig_I_pow = 0;          % slopes of I(q) in the selected q-region
            four_dat.fmt.fig_peak_pos = 0;       % peak position and amplitude
            % define positions
            screen_size = get(0,'screensize');
            %         four_dat.fmt.fig_pos_sym_int =  [150 150 screen_size(3)-300 screen_size(4)-300 ];
            %         four_dat.fmt.fig_pos_comb_all =  [150 150 screen_size(3)-300 screen_size(4)-300 ];
            % turn title off
            four_dat.fmt.title_on = 0;
            
            
            % if SAXS data has been corrected
            if plot_corrected_data
                processing.which_qrange_to_plot = 1;
                % calculate the fourier coefficients
                [four_dat] = calculate_four_dat(four_dat, one_projection.data);
            end
            % FOR SHALE (not in use)
            %         four_dat.ind_max_x = 71;
            %         four_dat.ind_max_y = 28;
            %         one_projection.par.ind_max_x = 71;
            %         one_projection.par.ind_max_y = 28;
            
            
            % if needed to apply a mask
            if (apply_mask)
                if ~(plot_corrected_data)
                    % adjust sizes, if mask and original data have different
                    % sizes
                    temp = [];
                    names = fieldnames(four_dat);
                    for pp = 2:9 % only the fields needed here
                        % y dimension
                        while size(four_dat.(names{pp}), 1) < size(one_projection.window_mask, 1)
                            temp = [four_dat.(names{pp})(1, :, :); four_dat.(names{pp})];
                            four_dat.(names{pp}) = temp;
                        end
                        % x dimension
                        while size(four_dat.(names{pp}), 2) < size(one_projection.window_mask, 2)
                            temp = [four_dat.(names{pp})(:, 1, :), four_dat.(names{pp})];
                            four_dat.(names{pp}) = temp;
                        end
                    end
                end
                % compile return structure for plotting and mask
                four_dat.f1_amp = four_dat.f1_amp.*one_projection.window_mask;
                four_dat.f2_amp = four_dat.f2_amp.*one_projection.window_mask;
                four_dat.f2_phase = four_dat.f2_phase.*one_projection.window_mask;
                four_dat.orientation = four_dat.orientation.*one_projection.window_mask;
                four_dat.degree_orient = four_dat.degree_orient.*one_projection.window_mask;
            end
            
            
            % plot the mesured symmetric intensity
            figh = figure(100);
            clf(figh)
            axesh = gca;
            subplot(2,3,1)
            imgh = imagesc(1);
            four_dat.fmt.figh = figh;
            four_dat.fmt.imgh = imgh;
            % plot
            plot_fourier_result(four_dat, processing);
            title(sprintf('Measured, rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y), 'FontSize', 10);
            
            
            % plot the measured orientation
            four_dat.fmt.fig_sym_int = 0;        % symmetric intensity over all segments: f1_amp
            four_dat.fmt.fig_asym_sym_int = 1;   % combination of orientation (colorwheel), degree of orientation (hue)
            figh = figure(100);
            subplot(2,3,4)
            imgh = imagesc(1);
            axesh = gca;
            four_dat.fmt.figh = figh;
            four_dat.fmt.imgh = imgh;
            % plot
            plot_fourier_result(four_dat, processing);
            title(axesh, sprintf('Measured, rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y), 'FontSize', 10);
            
            
            % Calculate 2D projections from the optimized data with SH
            filename = fullfile(base_path, sprintf('analysis/SASTT/%s/SH/%s/results/', sample_name, add_name));
            filename = fullfile(filename, sprintf('result_%s_q%d-%d_all_again_%s_reg_coeff_%0.3f.mat', sample_name, projection(1).par.r_sum{1}(1), ...
                projection(1).par.r_sum{1}(end), add_name,reg_coeff));
            load(filename)
            
            %%%%%%%%%%%%%%%%%%
            % For more efficiency the loop over projections should be done
            % inside SAXS_tomo_3D_err_metric
            [~,~,proj_SH] = SAXS_tomo_3D_err_metric([],p, s, one_projection);
            %%%%%%%%%%%%%%%%%%
            
            % load data as the variables must be reprocessed
            processing.load_data = 1; % 1 to force re-load and re-process
            processing.plot_data = 1;
            % define which q-range to plot
            processing.which_qrange_to_plot = 1;
            
            % calculate the fourier coefficients
            [four_dat] = calculate_four_dat(four_dat, proj_SH.projection);
            
            
            if (apply_mask)
                % compile return structure for plotting and mask
                four_dat.f1_amp = four_dat.f1_amp.*one_projection.window_mask;
                four_dat.f2_amp = four_dat.f2_amp.*one_projection.window_mask;
                four_dat.f2_phase = four_dat.f2_phase.*one_projection.window_mask;
                four_dat.orientation = four_dat.orientation.*one_projection.window_mask;
                four_dat.degree_orient = four_dat.degree_orient.*one_projection.window_mask;
            end
            
            
            % plot the simulated symmetric intensity
            four_dat.fmt.fig_sym_int = 1;        % symmetric intensity over all segments: f1_amp
            four_dat.fmt.fig_asym_sym_int = 0;   % combination of orientation (colorwheel), degree of orientation (hue)
            figh = figure(100);
            subplot(2,3,2)
            imgh = imagesc(1);
            four_dat.fmt.figh = figh;
            four_dat.fmt.imgh = imgh;
            % plot
            plot_fourier_result(four_dat, processing);
            title(sprintf('SASTT Sym. Int., rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y) , 'FontSize', 10);
            
            % plot the simulated orientation
            four_dat.fmt.fig_sym_int = 0;        % symmetric intensity over all segments: f1_amp
            four_dat.fmt.fig_asym_sym_int = 1;   % combination of orientation (colorwheel), degree of orientation (hue)
            figh = figure(100);
            subplot(2,3,5)
            imgh = imagesc(1);
            axesh = gca;
            
            four_dat.fmt.figh = figh;
            four_dat.fmt.imgh = imgh;
            % plot
            plot_fourier_result(four_dat, processing);
            title(axesh, sprintf('SASTT Orient., rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y), 'FontSize', 10);
            
            % include IRTT data
            filename = fullfile(base_path, sprintf('analysis/SASTT/%s/IRTT/%s/results/', sample_name, add_name));
            filename = fullfile(filename, sprintf('IRTT_%s_q_%d-%d_%s.mat', sample_name, projection(1).par.r_sum{1}(1), ...
                projection(1).par.r_sum{1}(end), add_name));
            load(filename)
            
            B = squeeze(B_segs(which_projections(kk),:,:));
            % projsize is [segments Y X]
            proj_IR = proj_simul_tensor(tomotensor, one_projection.Rot_exp, B);
            % change dimensions of proj_IR to [Y X segments]
            proj_IR = permute(proj_IR, [2 3 1]);
            % calculate  the fourier data
            [four_dat] = calculate_four_dat(four_dat, proj_IR);
            
            if (apply_mask)
                % compile return structure for plotting and mask
                four_dat.f1_amp = four_dat.f1_amp.*one_projection.window_mask;
                four_dat.f2_amp = four_dat.f2_amp.*one_projection.window_mask;
                four_dat.f2_phase = four_dat.f2_phase.*one_projection.window_mask;
                four_dat.orientation = four_dat.orientation.*one_projection.window_mask;
                four_dat.degree_orient = four_dat.degree_orient.*one_projection.window_mask;
            end
            
            % load data as the variables must be reprocessed
            processing.load_data = 1; % 1 to force re-load and re-process
            processing.plot_data = 1;
            % define which q-range to plot
            processing.which_qrange_to_plot = 1;
            
            % plot the simulated symmetric intensity
            four_dat.fmt.fig_sym_int = 1;        % symmetric intensity over all segments: f1_amp
            four_dat.fmt.fig_asym_sym_int = 0;   % combination of orientation (colorwheel), degree of orientation (hue)
            figh = figure(100);
            subplot(2,3,3)
            imgh = imagesc(1);
            four_dat.fmt.figh = figh;
            four_dat.fmt.imgh = imgh;
            % plot
            plot_fourier_result(four_dat, processing);
            title(sprintf('IRTT Sym. Int., rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y) , 'FontSize', 10);
            
            % plot the simulated orientation
            four_dat.fmt.fig_sym_int = 0;        % symmetric intensity over all segments: f1_amp
            four_dat.fmt.fig_asym_sym_int = 1;   % combination of orientation (colorwheel), degree of orientation (hue)
            figh = figure(100);
            subplot(2,3,6)
            imgh = imagesc(1);
            axesh = gca;
            
            four_dat.fmt.figh = figh;
            four_dat.fmt.imgh = imgh;
            % plot
            plot_fourier_result(four_dat, processing);
            title(axesh, sprintf('IRTT Orient., rotX %.1f, rotY %.1f', one_projection.rot_x, one_projection.rot_y), 'FontSize', 10);
        end
        
        % for saving the results
        if save_figure
            % save the figure
            %         figname = sprintf('scan_%05d_rotx%.0f_roty%.0f', four_dat.fnames.first_scan_no, one_projection.rot_x, one_projection.rot_y);
            figname = sprintf('rotx%.0f_roty%.0f', one_projection.rot_x, one_projection.rot_y);
            %make a new folder to save  the plotting comparison
            basedir_figures = fullfile(base_path,  sprintf('analysis/SASTT/%s/SH/%s/figures/%s', sample_name, add_name, figname));
            print(gcf, basedir_figures,'-dpng','-r300');
        end
        if save_video
            frame_plot = getframe(gcf);
            writeVideo(v, frame_plot);
            
        end
        pause(1)
        
        
        projSH_error(kk) = proj_SH.error;
        
        
        set(figh,'Position',[2026 155 1636 900]); %FKM
        
        
    end
    if save_video
        close(v)
    end
    
    
    savename_error = fullfile(base_path, sprintf('analysis/SASTT/%s/SH/%s/results/', sample_name, add_name));
    savename_error = [savename_error sprintf('projection_errors_reg_coeff_%0.3f.mat',reg_coeff)];
    save(savename_error,'projSH_error');
    
    %*-------------------------------------------------------------------------------------*
%|                                                                                     |
%|  Except where otherwise noted, this work is licensed under a                        |     
%|                                                                                     |   
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0                          |
%|  International (CC BY-NC-SA 4.0) license.                                           |
%|                                                                                     |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)                  |
%|                                                                                     |
%|      Author: CXS group, PSI                                                         |
%*------------------------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If this code, or subfunctions or parts of it, is used for research in a
%   publication or if it is fully or partially rewritten for another
%   computing language the authors and institution should be acknowledged
%   in written form and additionally you should cite:
%     M. Liebi, M. Georgiadis, A. Menzel, P. Schneider, J. Kohlbrecher,
%     O. Bunk, and M. Guizar-Sicairos, “Nanostructure surveys of
%     macroscopic specimens by small-angle scattering tensor tomography,”
%     Nature 527, 349-352 (2015).   (doi:10.1038/nature16056)
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Generate point plots

% point = [30 30];
% projNo = which_projections;
% 
% phi = linspace(180/32,180-180/32,16);
% figure(200);
% plot(phi,squeeze(projection(projNo).data(point(2),point(1),:)));
% hold on;
% plot(phi,squeeze(proj_SH.projection(point(2),point(1),:)));
% hold off;
% xticks([0 30 60 90 120 150 180]);
% ylabel('{\it I} (arb. units)');
% xlabel('\phi (deg.)');
% legend('Measured scattering','Simulated scattering','location','south');
% title(sprintf('Projection no. %d, point [%d, %d]',projNo,point(2),point(1)));
% axis square;
% ylim([0 inf]);
% xlim([0 180]);
% set(gcf,'color','w');
% 
% 
% %% Save plots for Figure S4.3 in paper 2 SI.
% figure(100);
% print(gcf,sprintf('/home/fredrimu/paper_figures/meas_sim_orientation_projNo_%d.tif',which_projections(1)),'-dtiff','-r600');
% 
% figure(200);
% print(gcf,sprintf('/home/fredrimu/paper_figures/I_phi_meas_sim_projNo_%d.tif',which_projections(1)),'-dtiff','-r600');