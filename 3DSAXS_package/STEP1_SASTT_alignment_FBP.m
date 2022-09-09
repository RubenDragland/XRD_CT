% This step aligns the projections. A tomogram is formed by N
% projections with tilt_angle = 0 based on filter back projection. The
% projections with different tilt angles (tilt_angle ~= 0) are aligned to 
% a sinthetic projection of the tomogram along two rotation angles. Indexes
% dx and dy are saved to align the projections during the optimization.


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

%error('Please run this script by section');

%% Step 1.1 load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_path = '~/Data10/';                    % = '~/Data10/' for online analysis,
                                            % provide the path for offline analysis
                                            % Ex: '/das/work/p16/p16649/'
sample_name = 'sample_name';                % name given in the saxs_caller_template
show_image = 1;
use_diode = 1; % 1: use diode (transmission data for alignment), 0 = use symmetric intensity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath utilities/
addpath utilities/tomo_alignment/
addpath(fullfile(base_path, 'cxs_software','base'));
addpath(fullfile(base_path, 'cxs_software','scanning_SAXS'));
addpath(fullfile(base_path, 'cxs_software', '3DSAXS'));
addpath(fullfile(base_path, 'cxs_software', 'tomo'));
addpath(fullfile(base_path, 'cxs_software', '3DSAXS', 'utilities'));
addpath(fullfile(base_path, 'cxs_software', '3DSAXS', 'utilities', 'tomo_alignment'));
file_name = sprintf('%sanalysis/SASTT/%s/projection_data/SASTT_%s.mat', base_path, sample_name, sample_name);
% load the projections
projection = load_projections(file_name);

% load the general parameters and paths
par = projection(1).par;
fnames = projection(1).fnames;
rot_axis = [];
tilt_axis = [];
%define the rotation axis according to the tomographic axis
for ii = 1:numel(projection)
    if par.tomo_axis_x == 0
        tilt_axis(ii) = projection(ii).rot_x;
        rot_axis(ii) = projection(ii).rot_y;
    else
    % Note, this option has not been tested
        tilt_axis(ii) = projection(ii).rot_y;
        rot_axis(ii) = projection(ii).rot_x;
    end
end


% select projections with tilt_axis = 0
index = find([tilt_axis]== 0);
%index = index([1:33]);  %% Remove projections (if needed) !write projections to keep

%stack measurements
stack_data = zeros(size(projection(1).diode, 1), size(projection(1).diode, 2), length(index));
for ii = 1:length(index)
    index_1 = index(ii);
    tilt_axis1 = tilt_axis(index_1);
    rot_axis1 = rot_axis(index_1);
    
    if use_diode
        data_select = -log10(projection(index_1).diode+min(projection(index_1).diode(projection(index_1).diode>0)));
        tit_name = 'Absorption';
    else
        %this has to be normalized by the maximum to work
        data_select = mean(projection(index_1).data, 3)./max(max(mean([projection(index).data], 3)));
        tit_name = 'Symmetric SAXS';
    end
    stack_data(:,:,ii) = data_select;
    % show images while processing the data
    if (show_image)
        figure(1)
        subplot(2, 4, [1, 2, 5, 6])
        imagesc(data_select)
        axis xy equal tight
        colormap bone
        title(sprintf('%s \nprojection %d/%d: tilt %d°, rot %d°',tit_name, ii, length(index),...
            tilt_axis1, rot_axis1));
        colorbar;
        caxis([0,1]);
        
        subplot(2, 4, [3,4])
        hold on
        errorbar(ii, mean(mean(data_select)), std(std(data_select)), 'ob');
        xlabel('Number of projections');
        ylabel('Mean pixel intensity')
        xlim([1, length(index)])
        grid on
        hold off
        
        subplot(2, 4, 7)
        hold on
        plot(ii, rot_axis1, 'ob');
        xlabel('Number of projections');
        ylabel('Rotation angles [°]')
        grid on
        xlim([1, length(index)])
        hold off
        
        if ii > 1
            subplot(2, 4, 8)
            hold on
            plot(ii, (rot_axis1 - rot_axis(index(ii-1))), 'ob');
            xlabel('Number of projections');
            ylabel('Angular step [°]')
            grid on
            xlim([1, length(index)])
            hold off
        end
        drawnow
    end
end


%% Step 1.2 Alignment in y (sample height)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp = 2;               % = 0 no display, =1 only final display, >1 every iteration
pixtol = 0.1;             % Tolerance of registration in pixels
rembias = true;         % Remove bias for y registration
maxorder = 1;           % Max order of bias to remove
limsy = [5 75];        % y-range for aligment. Assure that a common feature is within this range for all projections. 
deltaxal = 4;           % values along x. If too small the ROI becomes larger than the projection size. Best = 5
p.expshift = true; % Shift in phasor space at the end
params.interpmeth = 'linear';  % 'sinc' or 'linear' better for noise
params.alignx = false;  % Align horizontally with center of mass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

limsx=[1+deltaxal size(stack_data,2)-deltaxal];%THIS HAS TO BE CENTERED!!!
deltastack = ones(2, length(index))*0.1;

%make the tomographic axis vertical
% we might need to change for -y configurations
if par.tomo_axis_x
    stack_data1 = permute(stack_data, [2,1,3]);
else
    stack_data1 = stack_data;
end
[registration, ~] =  tomo.alignprojections_v4(stack_data1,limsy,limsx,deltastack,pixtol,rembias,maxorder,disp,params);

%% Step 1.3 refine vertical alignment excluding higher orders
maxorder = 2;
[registration, regstack1] =  tomo.alignprojections_v4(stack_data1,limsy,limsx,registration,pixtol,rembias,maxorder,disp,params);
% refine excluding higher orders
maxorder = 4;
[registration, regstack1] =  tomo.alignprojections_v4(regstack1,limsy,limsx,registration,pixtol,rembias,maxorder,disp,params);

%% Step 1.4 Horizontal alignment with tomographic consistency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pixtol = 0.01;              % pixel tolerance
disp = 2;                   % display the figures
paramsalign.interpmeth = 'linear';  % 'sinc' or 'linear'
paramsalign.filtertomo = 0.9;       % frequency cutoff
paramsalign.cliplow  = [ ];         % Clip air threshold
paramsalign.cliphigh = [-1e-3];     % Clip on sample threshold
paramsalign.binning = 0;            % Binning
paramsalign.usecircle = true;
initial_guess_offset = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

slice_num = round(size(regstack1, 1)/2);
sino = -squeeze(regstack1(slice_num,:,:)); % Create the sinogram
sino = real(utils.shiftpp2(sino,0.5,0)-utils.shiftpp2(sino,-0.5,0));
deltaslice = zeros(1,size(sino,2))+initial_guess_offset; % Create the sinogram;
tomo_axis = -([rot_axis(index)]);
[deltaslice, ~] =  tomo.alignslice_filt_v2(sino,tomo_axis,deltaslice,pixtol,disp,paramsalign);

% refinement of x alignemnt
paramsalign.cliphigh =  [-0.002];
paramsalign.filtertomo = 1;
[deltaslice, ~] =  tomo.alignslice_filt_v2(sino,tomo_axis,deltaslice,pixtol,disp,paramsalign);

%% Step 1.5 align multiple slices (10 slices considering the center of the projection)

slices = [round(size(regstack1, 1)/2)-5:round(size(regstack1, 1)/2)+5];
deltaxrefine = zeros(length(slices), length(index));
for ii = 1:length(slices)
    slices1 = slices(ii);
    display(['Aligning slice ' num2str(slices1)])
    sino = -squeeze(regstack1(slices1,:,:)); % Create the sinogram
    sino = real(utils.shiftpp2(sino,0.5,0)-utils.shiftpp2(sino,-0.5,0));
    [deltaaux, alignedaux] =  tomo.alignslice_filt_v2(sino,tomo_axis,deltaslice,pixtol,disp,paramsalign);
    deltaxrefine(ii,:) = deltaaux;
end

deltaxrefineav = mean(deltaxrefine,1);

fig_out = figure;
subplot(1,2,1)
imagesc(deltaxrefine);
xlabel(['Projection number'])
ylabel(['Slice number'])
title(['Displacements in x'])
%caxis([-5 5])
subplot(1,2,2)
plot(index, deltaxrefineav)
title(['Average displacements in x'])
xlabel(['Projection number'])

%% Step 1.6 shift the projections
interpmeth = 'linear';
regstack2 = regstack1*0;
for ii = 1:length(index)
    switch interpmeth
        case 'sinc'
            regstack2(:,:,ii) = angle(utils.shiftpp2(exp(1i*regstack1(:,:,ii)),0,deltaxrefineav(1,ii)));
        case 'linear'
            regstack2(:,:,ii) = angle(utils.shiftwrapbilinear(exp(1i*regstack1(:,:,ii)),0,deltaxrefineav(1,ii)));
    end
end
%% Step 1.7 see shifthed images
figure
for ii = 1:length(index)
    imagesc(regstack2(:,:,ii));
    title(['Projection ' num2str(ii) '/' num2str(length(index))]),
    colormap bone;
    axis xy equal tight;
    caxis([0 1])
    drawnow
    pause(0.1)
end
%% Step 1.8 make tomogram only for tilt_angle = 0
tomogram = zeros(size(regstack2, 2), size(regstack2, 2), size(regstack2, 1));
for i = 1:size(regstack2, 1) %along y
    sinogram = squeeze(regstack2(i,:,:)); % Create the sinogram
    freq_scale =1;
    data_tomogram_slice = tomo.iradonfast_v2(sinogram, tomo_axis, 'linear', 'Ram-Lak', ...
        size(regstack2,2), freq_scale);   % Calculate slice
    tomogram(:,:, i) = data_tomogram_slice;
    % sinogram
    figure(10);
    subplot(1,2,1)
    hold on
    imagesc(real(sinogram)')
    line([(size(sinogram,1)/2), (size(sinogram, 1)/2)] ,...
        [1, size(sinogram, 2)], 'Color','red','LineStyle','--')
    axis xy equal tight
    colormap bone
    colorbar
    title('Sinogram');
    xlabel(sprintf('Line integral'))
    ylabel('Projection''s angle \theta')
    
    % slice of the tomogram
    subplot(1,2,2)
    imagesc(data_tomogram_slice)
    hold on
    line([(size(data_tomogram_slice,2)/2), (size(data_tomogram_slice, 2)/2)] ,...
        [1, size(data_tomogram_slice, 1)], 'Color','red','LineStyle','--')
    line([1, size(data_tomogram_slice, 1)], [(size(data_tomogram_slice,2)/2), ...
        (size(data_tomogram_slice, 2)/2)],  'Color','red','LineStyle','--')
    hold off
    axis xy equal tight
    colormap bone
    axis off
    %caxis([0, 0.3])
    %colorbar
    title(sprintf('Slice %d\n along tomographic axis ',  i))
    drawnow
end

% correct  the tomogram for alignment problems
if par.tomo_axis_x
    % brings full tomogram [Z Y X] into beamline coordinate [X Y Z]
    % NOTICE THAT IT IS [Y X Z] in matlab [row, column, height]
    tomogram = permute(tomogram,[2,3,1]);
else
    % brings full tomogram [Z X Y] into beamline coordinate [X Y Z]
    % NOTICE THAT IT IS [Y X Z] in matlab [row, column, height]
    tomogram = permute(tomogram,[3,2,1]);
end

%% Step 1.9 Alignment of other tilt angles based on a projection of the tomogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cutting projections to avoid needles etc
fromx = [];
tox = [];
fromy = [];
toy = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

counter = 1;
for ii = 1:length(projection)
    fprintf('\n tilt = %d, rotation = %d/%d \n' , tilt_axis(ii),...
        rot_axis(ii), length(projection));
    % Size of tomogram
    N = size(tomogram);
    % defines the center in each direction
    % X appears in the second position of the meshgrid,
    % because meshgrid is (x,y,z) and size is (y,x,z)
    x = [1:N(2)] - ceil(N(2)/2); %
    y = [1:N(1)] - ceil(N(1)/2);
    z = [1:N(3)] - ceil(N(3)/2);
    
    % creates a mesh for arb_projection function
    [X, Y, Z] = meshgrid(x,y,z);
    
    %load the projections used for alignment
    if use_diode
        tosee = -log10(projection(ii).diode+eps);
    else % in case the sym_int is used
        tosee = mean(projection(ii).data, 3)./max(max(mean([projection.data], 3)));
    end
    Rot_exp_now = projection(ii).Rot_exp;
    
    % rotate the tomogram from step 2 based on the rotation matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EDIT HERE:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p.volume_upsampling = 1;  % Upsamples the volume by a factor of 2 to reduce artifacts
    p.method =  'bilinear';    % 'nearest' or 'bilinear' method for volume allocation
    p.filter_2D = 3;          % Strength of filter applied to image [0-3]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Register the measured projection and the projection from the tomogram
    xout = [1:size(tosee,2)]- ceil(size(tosee,2)/2);
    yout = [1:size(tosee,1)]- ceil(size(tosee,1)/2);
    % calculate the reference image based on the rotation matrix Rot_exp_now
    [proj_out_all, xout, yout] = arb_projection(tomogram,X,Y,Z,...
        Rot_exp_now,p,xout,yout);
    
    % reference image 
    ref_image =  proj_out_all(:,:);
    
    % define the images based on the output of diode or sym_int data
    to_be_aligned = tosee;
    
    % align the measured projection with the synthetic one
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EDIT HERE:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Fraction of a pixel for the fine search, default = 1
    % if > 1, then subpixel alignment takes place
    upsamp = 100;
    % displ     = 0  no information displayed (default)
    %           = 1  to display text information
    %           > 1  also shows images in figure(display)
    displ = 0;
    % Wfilt:     Fourier domain filter for registration (= 1 for no effect)
    Wfilt = 0.1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % choose cutout of image which corresond to each other, the differences have to
    % be added to the registration values (delta_x and delta_y)
    %[subim1, subim2, delta, deltafine, regionsout] = registersubimages_2...
    %(img1, img2, x1, y1, x2, y2, upsamp, displ,Wfilt)
    % choose the ROI that you want to align, since at high tilt angles the
    % full tomogram is not within the field of view
    % important in the case of needles appearing in the ROI
    % check some of the projections and adjust the ROI
    % if empty, it takes full size of projection
    
    [subim1, subim2, delta, ~,~] = utils.registersubimages_2...
        (ref_image, to_be_aligned, 1:length(xout), fromy:toy, ...
        fromx:tox, fromy:toy, upsamp, displ, Wfilt);
    subim2 = real(subim2);
    
    figure(354135)
    subplot(2,3,1)
    imagesc(subim1)
    axis xy equal tight
    title(sprintf('FBP projection \n tilt = %.2f, rotation = %.2f',...
        tilt_axis(ii), rot_axis(ii)));
    colormap jet
    %caxis([0, 1])
    colorbar
    
    subplot(2,3,2)
    imagesc(subim2)
    title(sprintf('Measurement \n tilt = %.2f, rotation = %.2f',...
        tilt_axis(ii), rot_axis(ii)));
    axis xy equal tight
    colormap jet
    % caxis([0, 1])
    colorbar
    
    % save results into structure
    projection(ii).dx = - double(delta(2));
    projection(ii).dy = - double(delta(1));
    
    %show difference between aligned and reference image
    subplot(2,3,3)
    rms_error = sqrt(((subim1-subim2).^2)./numel(subim1));
    imagesc(rms_error);
    colorbar
    axis xy equal tight
    title('RMS error');
    colormap jet
    
    subplot(2,3,[4:6])
    hold on
    errorbar(counter, mean(mean(rms_error)), std(std(rms_error)), 'ok')
    grid on
    xlabel('Number of projection')
    ylabel('Mean RMS error')
    xlim([0, length(projection)]);
    box on
    drawnow
    counter = counter +1;
end

%% Step 1.10 see aligned projections
figure
for ii = 1:length(projection)
    data = utils.shiftwrapbilinear(mean(projection(ii).data,3), - projection(ii).dy, - projection(ii).dx);
    subplot(1,2, 1) 
    imagesc(data);
    axis xy tight equal
    colormap bone
    title(sprintf('%d', ii))
     %caxis([0, 0.3]);
     subplot(1,2, 2)
     plot(ii, projection(ii).dx, 'or')
     hold on
     plot(ii, projection(ii).dy, 'ob')
     axis tight
     grid on
     legend([{'dx'}, {'dy'}])
     xlim([1, length(projection)])
      pause(0.1);
    drawnow;
end

%% Step 1.11 Remove needle and window mask
% Find not valid mask for needle or edge scattering from capillary
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

%%  create the masks with valid regions for 2D projections
% In the current implementation we remove three pixels from the projection,
% when the ASTRA toolbox is implemented this may not be needed, currently
% the projection on the last 2 pixels of the box is not reaiable.
mask_horiz = 2;                 % determine how much the FOV is going to shrink horizontally and vertically - Marios
mask_vert  = 2;                 % set = 0 to keep the whole FOV
threshold_valid = 20;            % After projecting the full 3D array with uniform = 1 value, whatever is below this threshold is considered not valid in the projection
[projection] = window_mask(projection, threshold_valid, mask_horiz, mask_vert);


%% Step 1.12 Show Figure: Check if the threshold for window and needle removal are correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stepping = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:stepping:length(projection)
    % plot the measurement window
    figure(3)
    clf;
    hold on
    subplot(1,2,1)
    imagesc(sum(projection(ii).window_mask,3))
    axis xy equal tight
    caxis([0,1]);
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
%% Step 1.13 saving the results
filename = fullfile(base_path,sprintf('analysis/SASTT/%s/projection_data/SASTT_%s_aligned_FBP.mat', ...
    projection(1).fnames.sub_dir, sample_name));
fprintf('Saving alignment to: %s \n', filename);
save(filename, 'projection', 'tomogram', '-v7.3');