% This step aligns the N projections with all rot and tilt angles simultaneously. A syntethic tomogram is
% reconstructed from the N projections. This is used to create a projection of the sinthetic tomogram along
% two rotation angles. Indexes dx and dy saved to align the projection in
% the optimization.


%*-------------------------------------------------------------------------------------*
%|                                                                                     |
%|  Except where otherwise noted, this work is licensed under a                        |
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

display('Please run this script by section');
return

addpath ..

%% STEP 1.1 Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent = cd ;
base_path = '/Data sets/'   %'~/Data10/';        % = '~/Data10/' for online analysis,
%parent = fullfile(parent, '..', '*' ) ;
base_path = [parent base_path] ;
                                % provide the path for offline analysis
                                % Ex: '/das/work/p16/p16649/'
sample_name = 'SASTT_carbon_knot_aligned_ASTRA_corrected'    %'sample_name';    % name given in the saxs_caller_template
sample_suffix='';               % suffix given to the data .mat file. Used for if the data was corrected, 
                                % for example if sgalil stages correction was introduced. 
n = [];                         % Apply some filtering to the diode data, n is the size of the filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the data % HAVE ALREADY ADDED NECESSARY DIRECTORIES TO MATLAB PATH
addpath utilities/
addpath utilities/tomo_alignment/
addpath(fullfile(base_path, 'cxs_software','base'));
addpath(fullfile(base_path, 'cxs_software','scanning_SAXS'));
addpath(fullfile(base_path, 'cxs_software', '3DSAXS'));
addpath(fullfile(base_path, 'cxs_software', '3DSAXS', 'utilities'));
addpath(fullfile(base_path, 'cxs_software', '3DSAXS', 'utilities', 'tomo_alignment'));
% file_name = sprintf('%sanalysis/SASTT/%s/projection_data/SASTT_%s%s.mat', base_path, sample_name, sample_name, sample_suffix); NOT MY FILE STRUCTURE
file_name = sprintf('%s%s.mat', base_path, sample_name); %CURRENTLY WHAT IS NEEDED.
% load the projections
projection = load_projections(file_name);

% create folder for saving the figures
% file_name = sprintf('%sanalysis/SASTT/%s/projection_data/figures/', base_path, sample_name);
file_name = sprintf('output%s', sample_name); % MAKE ONE DIR PER SAMPLE
if ~exist(file_name, 'file')
    mkdir(file_name);
end

if ~isempty(n)
    % This function multiplies by the transmission and then divides again
    % by a filtered transmission
    projection_filtered = filter_SASTT_diode(projection,n,'mean');
end


%% STEP 1.2 Remove projections (if needed): 
% if needed, run the lines below to check the projections
% imagesc(projection(21).diode); axis xy tight equal; colorbar
% imagesc(mean(projection(1).data, 3)); axis xy tight equal; colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
remove_proj = []; % = [], not to remove or = [56, 118] to remove projections 56 and 118, for example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correct projection
projection(1).par.remove_proj = []; %Initialize variable
if ~isempty(remove_proj)  && ~isfield(projection(1).par, 'remove_proj') 
    projection(1).par.remove_proj = remove_proj;
    projection(remove_proj) = [];
elseif (projection(1).par.remove_proj == remove_proj)
    fprintf('The projection %d has already been removed \n', remove_proj)
    return
end
%% STEP 1.2: compare the counts per projection, to see if air normalization is needed
diode_ct = [];
data_av = average_angular_sectors(projection);
saxs_ct = squeeze(mean(mean(data_av,1),2));
% select the number of projections
Nprojections = length(projection);

for ii = 1:Nprojections
   diode_ct = [diode_ct; (mean(projection(ii).diode(:)))];
end
figure(1)
clf
plot(1:Nprojections, saxs_ct, 'ko--')
hold on
xlabel('projection number')
ylabel('Mean SAXS Intensity [counts/s]')
yyaxis right

plot(1:Nprojections, diode_ct , '*--')
xlabel('projection number')
ylabel('Mean Diode Intensity [counts/s]')
legend({'mean saxs/10', 'mean diode'})
grid on
title('Compare the mean counts/projection: Is air normalization needed?')
hold off


%% STEP 1.3: normalize the SAXS signal by air signal per projection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norm_air = [];  
         % = [] : does not normalize
         % = 4 : number of air columns on the sides
plot_result = 1; % = 1 to plot result, otherwise = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% air normalization
if ~isempty(norm_air)
    [projection] = air_normalization(norm_air, projection, plot_result);
end

%% STEP 1.4 : Find not valid mask for needle or edge scattering from capillary
%  latter is only needed if something went wrong with the beam intensity 
%  during measurements for which bpm4i should be looked at but already the 
%  diode and the scattering as shown below could provide a clue
% It masks as non-valid parts of the projections that fall outside the reach 
% of the reconstructed volume, you can also include a threshold value to
% remove obscurations, such as sample mounting pin
% if needed, run the line below to check the threshold_needle value
% imagesc(projection(21).diode); axis xy tight equal; colorbar
% imagesc(mean(projection(1).data, 3)); axis xy tight equal; colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold_needle_max = [];    % usually around 0.9 for diode signal, = [] not to apply maximum
threshold_needle_min = [];    % use values from 0 to 1 for minimum threshold, or = [] not to apply minimum
p_needle_rem.use_diode_to_remove_needle = false;  % Either use diode data or the isotropic scattering
p_needle_rem.close_se_rad = 0;      % Radius of SE for closing operation, removes small dark spots from mask
p_needle_rem.erode_se_rad = 0;      % Radius of SE for erosion of mask, enlarges mask ( = [] to skip)
p_needle_rem.projections =  [];     % To which projections to apply the threshold ( = [] applies to all)
%%%% Needle bounds %%%%
% Some segmentations can be tricky, if this is the case you can use bounds
% to limit the threshold (before any morphological operation) to a subset
% of the data
%%% Sets to one above this height
p_needle_rem.max_needle_height = [];
% p_needle_rem.max_needle_height = size(projection(1).diode,1)*ones(numel(projection),1);
% p_needle_rem.max_needle_height(196:234) = 13;
%%% Sets to one to the right of this horizontal position
p_needle_rem.max_needle_right = [];
% p_needle_rem.max_needle_right = size(projection(1).diode,2)*ones(numel(projection),1);
% p_needle_rem.max_needle_right([207 244 246 247]) = 60;
%%% Sets to one to the left of this horizontal position
p_needle_rem.min_needle_left = [];
% p_needle_rem.min_needle_left = ones(numel(projection),1);
% p_needle_rem.min_needle_left(265:268) = 80;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[projection] = needle_removal(projection, threshold_needle_max, threshold_needle_min, p_needle_rem);


%% STEP 1.5 : Prepare/correct the data for alignment: either diode or SAXS for now
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_diode = 1;        % = 1: use diode (transmission data for alignment)
                      % = 0: use symmetric intensity: SAXS signal
show_movie = 1;       % = 1 to show  the images, = 0 not to              
stepping = 5; % define the step size to see less projections
frame_per_sec = 1;
axis_value = 'auto'; % scale the images, else 'auto'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This at some point should be removed and a more general description of 
% the sample orienation should be passed to the alignment code
data1 = [];
Nprojections = length(projection);
if use_diode
%      for ii = 1:Nprojections
%         data1(:,:,ii) = projection(ii).diode;
%      end
     data1 = pad_diode(projection);
     data1 = -log10(data1+min(data1(data1>0)));
else
    data1 = average_angular_sectors(projection);
end

% transposing for alignment
% Note: the signs of the angles for the case of tomo_axis_x remain to be tested
if projection(1).par.tomo_axis_x
    data1 = permute(data1,[2 1 3]);
end
for ii = 1:Nprojections
    if projection(1).par.tomo_axis_x
        par.tilt_angle(ii) = projection(ii).rot_y;
        theta(ii) =          projection(ii).rot_x;
    else
        par.tilt_angle(ii) = - projection(ii).rot_x;
        theta(ii) =          - projection(ii).rot_y;
    end
end

% add column: still needed for alignment
% [projection, data1] = add_column(projection, data1);


if isfield(projection(1),'needle_removal')&&(~isempty(projection(1).needle_removal))
    for ii = 1:Nprojections
        data(:,:,ii) = data1(:,:,ii).*projection(ii).needle_removal;
        if any(any(isnan(data(:,:,ii))))
            keyboard
        end
    end
else
    data = data1;
end

if show_movie
    figure(52)
    clf
    for ii = 1:stepping:Nprojections
        if isfield(projection(1),'needle_removal')&&(~isempty(projection(1).needle_removal))
            imagesc(data(:,:,ii).*projection(ii).needle_removal)
        else
            imagesc(data(:,:,ii))
        end
        axis equal tight xy
        colormap jet
        colorbar
        title(sprintf('Projection %d/%d \n Tilt %0.1f, Rot %0.1f', ii, Nprojections,...
            par.tilt_angle(ii), theta(ii)))
        colormap bone
        caxis(axis_value)
        pause(1/frame_per_sec)
    end
end

%% STEP 1.6 Alignment of projections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Npix = size(data, 2);  % X dimension
par.vert_range = [1:size(data, 1)]; % vertical layers for alignment. SAXS: whole y, diode: exclude needle
par.align_horizontal = true; % horizontal alignement
par.align_vertical = true; % vertical alignement: use only as second step is sample is difficult to align
par.high_pass_filter = 0.01; % Removes this fraction of lower spatial frequencies, recommended 0.01
par.min_step_size  = 0.001;  % how much subpixel precision is needed: 0.001
par.max_iter = 30; % maximal number of iterations: around 30
par.smooth_data = 0;  % smooth data with Gaussian filter (width given in pixels): used for weak scattering
par.show_diff_movie = false; % show difference after alignent
par.which_projection = [];  % to show specific projections, keep as [] to show all
par.stepping = 3;           % how many projections to skip if showing all, use 1 to show all
par.valid_angles = []; % Binary array for determining which projections are used in the reconstruction of the 
                            % tomogram used for alignment. The non-valid projections are not used for the tomogram but they are aligned, 
                            % useful when 45 degrees gives artifacts. Example: [1:20 50:size(data, 3)]; - do not use projections from 20 to 50.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create an input
shift = zeros(Nprojections,2);
%shift([21:49],2) = -10; Gives initial guess. Example here: shifts projections 21 to 49 by 10 pixel down. Used for example to align 45 degrees.
t_1 = tic;
for ii = 1:2
    % self consitency based alignment procedure based on the ASTRA toolbox
    [shift, par, tomogram1] = saxs.ASTRA_align_saxs(data, theta, Npix, shift, par);
    %par.align_vertical = true; % vertical alignement: in case it was false
    %for the first loop (ii = 1)
    % save image for comparison
    print(gcf, sprintf('%sanalysis/SASTT/%s/projection_data/figures/ASTRA_alignment_%s_%d', base_path, ...
        sample_name, sample_name, ii),'-dpng','-r500');
end
toc(t_1)
%% Additional refinement (optional)

% par.smooth_data = 0;  % smooth data with Gaussian filter (width given in pixels): used for weak scattering
%shift_aux = shift;
%shift_aux([21:49],2) = shift([22:49],2)-10; % In this example shifts projections 21 to 49 by 10 pixel down. Used for example to align 45 degrees.
%[shift, par, tomogram1] = saxs.ASTRA_align_saxs(data, theta, Npix, shift_aux, par);

%% STEP 1.7 Show the aligned projections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.show_diff_movie = true; % = 1 to show the differences
par.max_iter = 2;           % = number of iterations (minimum 2)
par.which_projection = [];  % to show specific projections, keep as [] to show all
par.stepping = 1;           % how many projections to skip if showing all, use 1 to show all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[shift, par, tomogram1] = saxs.ASTRA_align_saxs(data, theta, Npix, shift, par);
% save image for comparison
print(gcf, sprintf('%sanalysis/SASTT/%s/projection_data/figures/ASTRA_offset_%s', base_path, ...
    sample_name, sample_name),'-dpng','-r300');
%% STEP 1.8 Save the projection shifts
% permute the tomogram
tomogram = permute(tomogram1, [3, 1, 2]);
tomogram = flip(tomogram, 3);
% shift the aligned sample
for ii = 1:size(data,3)
    projection(ii).dx =   shift(ii, 1);
    projection(ii).dy =   shift(ii, 2) ;
    corr_data(:,:,ii) = utils.shiftwrapbilinear(data(:,:,ii), -projection(ii).dy, -projection(ii).dx);
end
% check the reconstruction/alignment
view3(tomogram)

%% STEP 1.9 Show the aligned sample and plot the displacements
figure(14564)
clf
subplot(1,2, 1)
    plotting.imagesc3D(corr_data(:,:,:));
    axis xy tight equal
    colormap bone
    title('Projection %d');
    %caxis([0, 0.3]);
    subplot(1,2, 2)
    plot([projection(:).dx], 'or')
    hold on
    plot([projection(:).dy], 'ob')
    axis tight
    grid on
    legend([{'dx'}, {'dy'}])
    xlim([1, Nprojections])
% save figure for comparison
print(gcf, sprintf('%sanalysis/SASTT/%s/projection_data/figures/ASTRA_dx_dy_%s', base_path, ...
    sample_name, sample_name),'-dpng','-r300');

%% STEP 1.10 create the masks with valid regions for 2D projections
% In the current implementation we remove three pixels from the projection,
% when the ASTRA toolbox is implemented this may not be needed, currently
% the projection on the last 2 pixels of the box is not reaiable.
mask_horiz = 2;                 % determine how much the FOV is going to shrink horizontally and vertically
mask_vert  = 2;                 % set = 0 to keep the whole FOV
threshold_valid = 20;            % After projecting the full 3D array with uniform = 1 value, whatever is below this threshold is considered not valid in the projection
[projection] = window_mask(projection, threshold_valid, mask_horiz, mask_vert);

%% STEP 1.11 Show mask: Check if the threshold for window and needle removal are correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stepping = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:stepping:length(projection)
    % plot the measurement window
    figure(3)
    clf;
    
    hold on
    subplot(1,2,1)
    imagesc(sum(projection(ii).window_mask,3))
    axis xy equal tight
    title(sprintf('Define the window \n %d/%d',...
        ii, length(projection)));
    
    % plot_visible projection
    subplot(1,2,2)
    if use_diode
        imagesc((sum(projection(ii).window_mask,3)).*projection(ii).diode)
    else
        imagesc((sum(projection(ii).window_mask,3)).*mean(projection(ii).data, 3))
    end
    axis xy equal tight
    colormap bone
    title(sprintf('Region considered \n  %d/%d',...
        ii, length(projection)));
    pause(0.1)
end

%% STEP 1.12 Save the aligned projections
filename = fullfile(base_path,sprintf('analysis/SASTT/%s/projection_data/SASTT_%s_aligned_ASTRA.mat', ...
    sample_name, sample_name));
fprintf('Saving alignment to: %s \n', filename);
save(filename, 'projection', '-v7.3','-nocompression');

%% EXTRA: Compare projections at (0,0) angles 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im_columns = 4; % in how many columns should the data be shown?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = find(par.tilt_angle == 0 & theta == 0);

% get min and max value to equalize caxis in the end                  
temp=[];
for ii = 1:length(index)
    temp{ii}=data(:,:,index(ii));
end
    max0=max([temp{:}]);
    min0=min([temp{:}]);     
    
% compare projections at tilt and rotation angle = 0
figure
for ii = 1:length(index)
    subplot(round(length(index)/im_columns), im_columns, ii)
    %imagesc(mean(projection(index(ii)).data, 3))
    imagesc(temp{ii})
    axis equal tight xy
    colorbar
    title(sprintf('Projection %.0f', index(ii)))
    axis off
    colormap jet
    %pause(1)
    caxis([min(min0) max(max0)])
end

% get min and max value to equalize caxis in the end    % Marios                     
temp=[];
for ii = 2:length(index)
    temp{ii}=data(:,:,index(ii))-data(:,:,index(ii-1)); 
end
    max0=max([temp{:}]);                                  % Marios
    min0=min([temp{:}]);     


% difference plot between the projections
figure
for ii = 2:length(index)
    subplot(round(length(index)/im_columns), im_columns, ii-1)
    %imagesc(mean(projection(index(ii)).data, 3))
    imagesc(temp{ii})
    axis equal tight xy
    colorbar
    title(sprintf('Projection %.0f-%.0f', index(ii), index(ii-1)))
    axis off
    colormap jet
    %pause(1)
    caxis([min(min0) max(max0)])
end
temp=[];

%% EXTRA: only needed in case the range in y should be reduced
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % EDIT:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y_range = [];                               % = [] optimizes the whole y dimension, 
%                                             % = [25:70], for example, to select a range in y
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                             
% % need to select a range in y
% if ~isempty(y_range)
%     for ii = 1:Nprojections
%         % for this sample
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         projection(ii).data = projection(ii).data(y_range, :, :);
%         projection(ii).diode = projection(ii).diode(y_range,:);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     end
% end

