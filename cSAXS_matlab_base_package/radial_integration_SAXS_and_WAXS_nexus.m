% radial_integration_SAXS_and_WAXS.m
% Template for radial integration made around 2015
% Changes:
% 2016-08-22: define mask files at the beginning, allowing for a flag in case it needs to be repeated  
% add the save fast and v6

% License at the end of script

clear all
close all

%% step 0: add the path for the matlab-scripts (fill in userID,detno and specdatfile)
addpath .

userID  = [beamline.identify_eaccount '_']; %e-account followed by underline

detector = 'pilatus_1';    % detector name: 'pilatus_1' for SAXS Pilatus 2M, 'pilatus_2' for WAXS Pilatus 300k,
                           % 'eiger' for SAXS Eiger 500 k
save_format = '-v6'; % which data format to save? '-v6' is the standard.
file_flag=''; % suffix for filenames for valid pixel mask, beamstop mask coordinates and integration mask
              % this suffix is also used for radial_integration folders 
              % Example: '_2M_at_two_meters' (include underscore recommended)
              % Leave empty '' for default folder and filenames.
homedir  = '~/Data10/';

prepare_mask = true;

% homedir  = '/das/work/p16/p16268/'; % change here for offline analysis

%%  Scan numbers of standards
%glassy carbon, glassy carbon moved detector to side, air scattering, first
%one is glassy carbon used to remove beamstop later
 scannr = [3191 3191 3190];
 %AgBE (for SAXS and WAXS), LaB6 (for WAXS), Si (for WAXS)
   todo      = [3203 3192 3193];
  legendstr = {'AgBE';'LaB6';'Si'};


%%%%%%%% Modify until here %%%%%%%%

if (strcmp(detector,'pilatus_1'))||(strcmp(detector,'eiger'))
    Ifsaxs = 1; % =1 for SAXS detector, =0 for WAXS
elseif strcmp(detector,'pilatus_2')
    Ifsaxs = 0;
else
    error('\nUnknown detector %s, please check.',detector);
end

addpath(fullfile(homedir,'cxs_software','base'));

%CHANGE: spec dat file
SpecDatFile = homedir; 

datadir = fullfile(homedir,'data');

if strcmp(detector,'pilatus_2')
    integdir = fullfile(homedir,'analysis', sprintf('radial_integration_waxs%s/',file_flag));
elseif strcmp(detector,'pilatus_1')
    integdir = fullfile(homedir,'analysis', sprintf('radial_integration%s/',file_flag));
elseif strcmp(detector,'eiger')
    integdir = fullfile(homedir, 'analysis', sprintf('radial_integration_eiger%s/',file_flag));
end
if strcmp(detector,'pilatus_2')
    outdir   = fullfile(homedir,'analysis',sprintf('data_waxs%s',file_flag));
elseif strcmp(detector,'pilatus_1')
    outdir   = fullfile(homedir,'analysis', sprintf('data%s', file_flag));
elseif strcmp(detector,'eiger')
    outdir   = fullfile(homedir,'analysis',sprintf('data_eiger%s',file_flag));
end

maskfilename = fullfile(outdir,sprintf('%s_valid_mask%s.h5', detector,file_flag)); 
integmaskfilename=fullfile(outdir,sprintf('%s_integration_masks%s.h5',detector,file_flag)); 

integrate_range_args = {'Detector', detector};

dirs = whos('-regexp','.*dir$');
for ii=1:numel(dirs)
  dir_to_do = eval(dirs(ii).name);
  if ~exist(dir_to_do,'dir')
    fprintf('creating directory %s\n', dir_to_do);
    system(sprintf('mkdir -p %s',dir_to_do));
  end
end


%% step 1: prepare the valid pixel mask

if prepare_mask
  fprintf('preparing the valid pixel mask\n');

% calculating the union of several valid pixel masks
%   starting with a rather dark file to discriminate hot pixels
  system(sprintf('rm -f %s', maskfilename));
 
  for ii=scannr
    beamline.prep_valid_mask_nexus(fullfile(homedir,'data',utils.compile_x12sa_dirname(ii)), ...
        detector, ...
        'ThresholdDark',1, ...
        'ThresholdHot',20, ...
        'Extend','or', ...
        'FilenameValidMask',maskfilename, ...
        'AdjustOrientation', false);
      %  'FigNo',ii==scannr(end));
  end  
end
%% step 2: cut out beam stop and shadows manually (for WAXS only necessary if there is a shadow)
% Instructions:
% 1) Run this section and patiently wait until two figures pop up. The last
%    one is called "Mask selection"
% 2) In "Mask selection" select the shape you want, either circle,
%    rectangle, or polygon.
% 3) Click on "Add more pixels".
% 4) Go to figure 555 where you can select the region you want to remove
%    from the valid mask, for polygon make sure to close the polygon.
% 5) Double click inside the selected region
% 6) You will see the region has been added on "Mask selection". Click on
%    "add more pixels" to add another region, or click "Quit" to exit and
%    save the mask.
% 7) The first time you create a mask it will save the history in the analysis/data
%    directory. IN this way, the next time that you make a mask you can
%    have access to the history by pressing the "load" button. THe history
%    also saves the layers separately.
% 8) If you press "reset" it will start from scratch a new mask and one can
%    load the latest mask which was created (typically the bad pixels
%    without the beamstop) or the last history that was created (which
%    includes different layers)

save_mask_history = true;
overwrite_history = true;

filename = fullfile(homedir,'data',utils.compile_x12sa_dirname(scannr(1)),sprintf('%s1_%05d.h5',userID,scannr(1)) );
frame = io.nexus_read(filename, 'filter', detector);

if save_mask_history
    [fpathHist, fnameHist] = fileparts(maskfilename);
    if overwrite_history
        histSuffix = '';
    else
        histSuffix = ['_' datestr(datetime(), 30)]; % add date and time to history file (ISO 8601: yyyymmddTHHMMSS)
    end
    fnameHistFull = fullfile(fpathHist, [fnameHist '_history' histSuffix '.h5']);
else
    fnameHistFull = [];
end

if prepare_mask
    figure(555);
    ax = imagesc(log10(abs(frame.data(:,:,1)))); axis equal tight xy
    colormap(plotting.colormaps.magma);
    colorbar();
    pause(3)
    beamstop_mask_plot = beamline.create_mask('axis', ax, 'disable_save', true, 'fname_history', fnameHistFull, 'orientation', [frame.orientation.transpose frame.orientation.rot90], 'load_mask', maskfilename);
    beamstop_mask_save = math.applyTransform(beamstop_mask_plot, ...
        [frame.orientation.transpose frame.orientation.rot90], 'transpose_rot90', true);
    figure(556)
    imagesc(log10(abs(beamstop_mask_plot.*frame.data(:,:,1)))); axis equal tight xy;
    colormap(plotting.colormaps.magma);
    colorbar();
    
    io.save_mask_nexus(maskfilename, beamstop_mask_save);
    fprintf('Saving mask in %s\n',maskfilename);
    
    close(555)
end

%% show silver behenate scattering to find the radius of the first ring (only SAXS)

filename = fullfile(homedir,'data',utils.compile_x12sa_dirname(todo(1)),sprintf('%s1_%05d.h5',userID,todo(1)));
[frame, valid_mask] = io.nexus_read(filename, 'filter', detector, 'mask', maskfilename);

frame_plot = log10(frame.data(:,:,1).*valid_mask);
fig_ag = figure(1);
imagesc(frame_plot); axis equal tight xy
colormap(plotting.colormaps.viridis);
colorbar();
valmax = math.sp_quantile(frame_plot, [0 1-1e-3],10);
caxis([-1 valmax(2)]);
title(sprintf('Silver behenate - Scan %05d', todo(1)));
set(gcf,'Name','Silver behenate');

%% here you have to give some manual inputs to run step 3
% for SAXS you have to put y pixel value of the the silver behenate ring above the beamstop, and the order of the peak that you chose
if Ifsaxs
    order_AgBE = 2;
    y_from = 1177;
    y_to   = 1182;
    cen_guess = [];  %[y,x] ; leave empty, i.e. cen_guess=[], for automatic guess;
    %and choose how many sectors you want to do the integration (16 for
    %anisotropic scattering, 1 for isotropic scattering
    num_segments=16;
    skip_center_refinement = false;
   
else
    %for WAXS you can run with the default values to start with and adjust in
    %case an error appears or the fit (shown in figure 4) is bad
    
    openfig('+beamline/WAXS_standards.fig','reuse');
    %give the order of the first silver behenate ring appearing
    %(compare with WAXS_standards.fig)
    order_AgBe=6;
    
    %parameter used in finding the x-position, default 5, if in figure 20 the
    %blue curve is all zeros, lower this value (necessary for low intensity of
    %silver behenate measurement
    d = 4;
    
    %threshold to find WAXS peak of standards, default is 50, might be lowered
    %for lower intensities
    threshold=[2 10 100];
    %if wrong peaks are found tune finding the right peaks with the window
    %where peaks are being searched here, default is min=0 and max=1500,
    %(see WAXS_standards.fig)
    min_AgBE=0;
    max_AgBE=480;
    min_Si=0;
    max_Si=1500;
    min_LaB6=500;
    max_LaB6=1400;
    
end
% step 3: prepare integration mask
% For the WAXS mask this is still a bit clunky. You can adjust above the
% min and max values where it will look for a peak and the threshold. Also
% in the fit for the horizonal position make sure there is both red and
% blue peaks for the fitting, if not you can adjust the d parameter above.
% Decreasing it helps when the silver behenate scattering is low.

filename = fullfile(homedir,'data',utils.compile_x12sa_dirname(todo(1)),sprintf('%s1_%05d.h5',userID,todo(1)));
nexus_args = {'H5Location',['/entry/instrument/' detector '/data']};
% Arguments needed to convert coordinates from beamline orientation to raw data orientation
convert_location_args = {'framesz', size(frame.data(:,:,1)), 'vecIn', ...
    [frame.orientation.transpose frame.orientation.rot90], 'vecTypeIn', 'transpose_rot90'};


S = io.nexus_read(homedir,'scanNr',todo(1),'filter','spec');

[I, mask] = io.nexus_read(filename, 'filter', detector, 'mask', maskfilename, 'orientation', [0 0]);
I = mean(I.data,3).*mask;

if Ifsaxs
    while true
        if isempty(cen_guess)
            J = ifftn(fftn(I,size(I)*2-[1 1]).^2);
            cen_guess_raw = math.peakfit2d(J)/2; %[y,x]
            cen_guess = math.convert_location('pos', cen_guess_raw, ...
            'framesz', size(I), 'vecIn', [0 0], ...
            'vecTypeIn', 'transpose_rot90', 'vecOut', ...
            [frame.orientation.transpose frame.orientation.rot90], ...
            'vecTypeOut', 'transpose_rot90');
        else
            cen_guess_raw = math.convert_location('pos', cen_guess, ...
                'framesz', size(frame.data(:,:,1)), 'vecIn', [frame.orientation.transpose frame.orientation.rot90], ...
                'vecTypeIn', 'transpose_rot90');
        end
        
        filename = fullfile(homedir,'data',utils.compile_x12sa_dirname(todo(1)),sprintf('%s1_%05d.h5',userID,todo(1)));
        
        %y_from-cen_guess_raw(1)
        coord_from_raw = math.convert_location('pos', [y_from cen_guess(2)], ...
                'framesz', size(frame.data(:,:,1)), 'vecIn', [frame.orientation.transpose frame.orientation.rot90], ...
                'vecTypeIn', 'transpose_rot90');
        radiusFrom = sqrt(sum((cen_guess_raw(:) - coord_from_raw(:)).^2));
        
        %y_to-cen_guess_raw(1)
        coord_to_raw = math.convert_location('pos', [y_to cen_guess(2)], ...
            'framesz', size(frame.data(:,:,1)), 'vecIn', [frame.orientation.transpose frame.orientation.rot90], ...
            'vecTypeIn', 'transpose_rot90');
        radiusTo = sqrt(sum((cen_guess_raw(:) - coord_to_raw(:)).^2));
        if radiusTo < radiusFrom
            error('radiusTo < radiusFrom, please change above y_to and y_from')
        end
        
        % Notice "cen" is in raw-frame-orientation coordinates
        if ~skip_center_refinement
            [cen]=utils.get_beam_center_nexus(filename,detector,'GuessX',cen_guess_raw(2),'GuessY',cen_guess_raw(1), ...
                'RadiusFrom',radiusFrom,'RadiusTo', radiusTo, ...
                'TestX',4,'TestY',4,'FilenameValidMask',maskfilename, nexus_args{:});
        else
            fprintf('Skipping center refinement\n')
            cen = flipud(cen_guess_raw(:)); %flipping needed because of different xy or yx conventions in different codes
        end
        
        cen_print = math.convert_location('pos', [cen(1), cen(2)], ...
            'framesz', size(frame.data(:,:,1)), 'vecIn', [0 0], ...
            'vecTypeIn', 'transpose_rot90', 'vecOut', ...
            [frame.orientation.transpose frame.orientation.rot90], ...
            'vecTypeOut', 'transpose_rot90');
        if ~strcmpi(input(sprintf('The center (Y/X) changed from (%.3f,%.3f) to (%.3f,%.3f).\nDo you want to repeat the search with the new guess? Y/N [Y]: ', cen_guess(1), cen_guess(2), cen_print(2), cen_print(1)), 's'), 'n')
            cen_guess = [cen_print(2) cen_print(1)];
        else
            break;
        end
    end
else
    % this isn't nice yet
    % i)   it depends on the chosen orientation on how to read
    %      detector-2 images
    % ii)  it merely finds maximum values instead of fitting, possibly
    %      with sub-pixel precision
    % iii) as a consequence, figuring out which values are trustworthy
    %      is done rather crudly
    
    %%% First determine center line parameter of diffraction patterns on the nexus
    %%% orientation, done like this since the code below is heavily
    %%% dependant on the pattern orientation
    [I_nexus_orient, mask_nexus_orient] = io.nexus_read(filename, 'filter', detector, 'mask', maskfilename);
    I_nexus_orient = mean(I_nexus_orient.data,3).*mask_nexus_orient;
    dx = 30;
    [s1,s2] = size(I_nexus_orient);
    J = ifft(fft(I_nexus_orient,s1*2-1,1).^2,[],1);
    [~,n] = max(J);
    w = std(I_nexus_orient,1,1)./sqrt(mean(I_nexus_orient,1));
    o = 1:numel(n);
    o = o(w>d);
    n = n(w>d)/2;
    o = o(abs(n-s1/2)<dx);
    n = n(abs(n-s1/2)<dx);
    x = s1/2+linspace(-dx,dx,4*dx+1);
    m = histc(n,x);
    [~,n0] = max(m);
    %%% fitting curve to find the center
    s = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[  0,s1/2-dx,  0,  0,  0],...
        'Upper',[Inf,s1/2+dx,Inf,Inf,Inf],...
        'Startpoint',[10,x(n0),1,10,1]);
    f = fittype('a*exp(-((x-b)/c)^2)+d*exp(-((x-n)/e)^2)', ...
        'problem','n','options',s);
    [c,~] = fit(x',m',f,'problem',s1/2);
    
    if all(m==0)
        error('Could not find signal to determine the center in x, please lower the value of d')
    end
    
    figure(20)
    plot(x,m)
    hold on
    plot(c,'r');
    hold off
    title('Finding center in x')
    cen1 = c.b;
    o = o(abs(n-cen1)<=1);
    n = n(abs(n-cen1)<=1);
    
% Pretty useless figure
%     figure(10)
%     hold on
%     plot(o,n,'w.')
%     plot([1 s2],[1 1]*round(cen1),'w')
%     x = 1:s2;
%     plot(x(mask_nexus_orient(round(cen1),:)>0), ...
%         log(I_nexus_orient(round(cen1),mask_nexus_orient(round(cen1),:)>0))/ ...
%         max(log(I_nexus_orient(round(cen1),mask_nexus_orient(round(cen1),:)>0)))*s1, ...
%         'k')
%     hold off
    
    figure(fig_ag)
    plotting.hline(cen1,'-c')
    
    figure(30)
    WAXS = zeros(s2,numel(todo));
    WAXS(:,1) = I_nexus_orient(round(cen1),:);
    for ii=2:numel(todo)
        I_nexus_orient = io.nexus_read(homedir,'scanNr',todo(ii),'filter',detector);
        I_nexus_orient.data(I_nexus_orient.data<0) = 0;
        % What about the valid mask? Currently not used here
        WAXS(:,ii) = mean(I_nexus_orient.data(round(cen1),:,:),3);
    end
    
    h = semilogy(WAXS);
    legend(legendstr)
    
    % finding peaks "automatically"
    x_coord = [];
    q_coord = [];
    hold on
    peaks = cell(1,size(WAXS,2));
    for ii=1:size(WAXS,2)
        %the treshhold value, default set to 50, might be adjusted
        peaks{ii} = utils.peakfinder((WAXS(:,ii)),threshold(ii));
        %peaks{ii} = peakfinder((WAXS(:,ii)),50);
        if strcmp(legendstr{ii},'AgBE')
            tmp = peaks{ii};
            tmp = tmp(tmp>=min_AgBE);
            peaks{ii} = tmp(tmp<=max_AgBE);
            
        end
        if strcmp(legendstr{ii},'Si')
            tmp = peaks{ii};
            tmp = tmp(tmp>=min_Si);
            peaks{ii} = tmp(tmp<=max_Si);
        end
        if strcmp(legendstr{ii},'LaB6')
            tmp = peaks{ii};
            tmp = tmp(tmp<=max_LaB6);
            peaks{ii} = tmp(tmp>=min_LaB6);
            
        end
        x_coord = vertcat(x_coord,peaks{ii});
        if strcmp(legendstr{ii},'AgBE')
            q0 = 2*pi/58.38;
            q_coord = horzcat(q_coord,q0*(order_AgBe+(0:numel(peaks{ii})-1)));
        elseif strcmp(legendstr{ii},'LaB6')
            q0 = 2*pi/4.1549;
            q_coord = horzcat(q_coord,q0*sqrt((1:numel(peaks{ii}))));
        elseif strcmp(legendstr{ii},'Si')
            q0 = 2*pi/5.4308;
            q_coord = horzcat(q_coord,q0*sqrt(3));
        end
        semilogy(peaks{ii},WAXS(peaks{ii},ii),'.', ...
            'Color',get(h(ii),'Color'), ...
            'MarkerSize',24)
    end
    hold off
    
    figure(40); clf
    if (numel(x_coord)>3)
        %       fprintf('%f\t%f\n',[x_coord';q_coord])
        %       % a
        %       % b
        %       % c
        s = fitoptions('Method','NonlinearLeastSquares',...
            'Lower'     ,[-Inf,-Inf,  0],...
            'Upper'     ,[ Inf,   0,1e3],...
            'Startpoint',[s2/2, 200,550]);
        f = fittype('4*pi/l*sin((atan((a-b)*p/c)+atan((x-a)*p/c))/2)', ...
            'problem',{'p','l'},'options',s);
        [c,~] = fit(x_coord,q_coord.',f,'problem',{.172,12.398/double(S.mokev)});
        subplot(2,1,1)
        plot(x_coord,q_coord,'x');
        ylabel('Tabulated values for peaks [A^{-1}]')
        xlabel('pixel value of found peaks')
        hold on
        drawnow;
        tmp = axis;
        x = linspace(c.b,tmp(2));
        plot(x,feval(c,x),'r');
        subplot(2,1,2)
        bar(x_coord,feval(c,x_coord)-q_coord');
        xlim(tmp(1:2));
        dc = confint(c); % Confidence intervals
        dc = (dc(2,:)-dc(1,:))/2;
        fprintf(['detector distance:\t%.1fmm,    \t%.1fmm\n', ...
            'center of rings:  \t%.1fpixels,\t%.1fpixels\n', ...
            'angle of detector:\t%.1fdeg,   \t%.1fdeg.\n'],   ...
            c.c,dc(3), ...
            c.b,dc(2), ...
            atan((c.a-c.b)*c.p/c.c)/pi*180, ...
            180/pi*c.p/c.c*sqrt(dc(1)^2+dc(2)^2 + ((c.a-c.b)/c.c*dc(3))^2));
    end
end


tic
filename = fullfile(homedir,'data',utils.compile_x12sa_dirname(todo(1)),sprintf('%s1_%05d.h5',userID,todo(1)));
S = io.nexus_read(filename, 'filter', 'spec');
wavelength = 12.3984197/S.mokev;  % In Angstrom
if Ifsaxs
    fprintf('preparing the integration mask(s)\n');
    
    % Notice "cen" is already in raw-frame-orientation coordinates
    beamline.prep_integ_masks_nexus(filename, cen, 'Det',detector, ...
        'Saxs', Ifsaxs, ...
        'NoOfSegments',num_segments, ...
        'FilenameValidMask',maskfilename, ...
        'FilenameIntegMasks',integmaskfilename, nexus_args{:},'ConvertLocationArgs',convert_location_args);
    
    beamline.integrate_range_nexus(todo(1),'stopScanNr',todo(1),...
        'detector', detector,'outdirData',integdir,'basePath',datadir,'integMask',integmaskfilename);
    
else
    
    %%% Convert [c.b cen1] and [c.a cen1] to raw-frame-orientation
    %%% coordinates
    cen_beam_waxs_raw_yx = math.convert_location('pos', [cen1 c.b], ...
        'framesz', size(frame.data(:,:,1)), 'vecIn', [frame.orientation.transpose frame.orientation.rot90], ...
        'vecTypeIn', 'transpose_rot90');

    normxy_beam_waxs_raw_yx = math.convert_location('pos', [cen1 c.a], ...
        'framesz', size(frame.data(:,:,1)), 'vecIn', [frame.orientation.transpose frame.orientation.rot90], ...
        'vecTypeIn', 'transpose_rot90');
    
    fprintf('preparing the integration mask(s)\n');
    %%% notice I gotta swap the order of the center because they are in yx
    %%% and this function receives xy
    beamline.prep_integ_masks_nexus(filename, [cen_beam_waxs_raw_yx(2) cen_beam_waxs_raw_yx(1)], 'Det',detector, 'Saxs', Ifsaxs, ...
        'Wavelength', wavelength, 'NormalXY', [normxy_beam_waxs_raw_yx(2) normxy_beam_waxs_raw_yx(1)], 'DetDist_mm', c.c, ...
        'PixelSize_mm', frame.x_pixel_size*1e-3, 'NoOfSegments',1, 'FilenameValidMask',maskfilename, ...
        'FilenameIntegMasks',integmaskfilename, 'DisplayValidMask',0,...
        nexus_args{:},'ConvertLocationArgs',convert_location_args);
end
toc


%% calculate detector distance (SAXS only) check in Figure 100 if the peak_agbe really is the 1st order AgBE
if Ifsaxs
    [x,y] = plotting.plot_radial_integ(sprintf('%s%s%d_%05d_00000_00000_integ.h5',integdir,userID,1,todo(1)));
    %%the 1st order silver behenate is at ... pixels
    %peakfinder(log(y(10:end)),1);
    peaks2 = utils.peakfinder(log(y(10:end)),1);
    peak_agbe = x(peaks2(order_AgBE+1))+9 %normally the 1st order AgBE, check!
    detector_distance = peak_agbe*.172/tan(2*asin(wavelength*order_AgBE/(2*58.38)))
    plotting.plot_radial_integ(sprintf('%s%s%d_%05d_00000_00000_integ.h5',integdir,userID,1,todo(1)),'FigNo',101,'SegAvg',0);
end
%% redo SAXS integration mask now it will take the detector distance into account and also save the q-value
if Ifsaxs
    if (strcmp(detector,'pilatus_1'))
        detector_pixelsize = 0.172;
    elseif (strcmp(detector,'eiger'))
        detector_pixelsize = 0.075;
    end
    filename = fullfile(homedir,'data',utils.compile_x12sa_dirname(todo(1)),sprintf('%s1_%05d.h5',userID,todo(1)));
    S = io.nexus_read(filename, 'filter', 'spec');
    fprintf('preparing the integration mask(s)\n');
    beamline.prep_integ_masks_nexus(filename, cen, 'Det',detector, 'Saxs',Ifsaxs,...
        'NoOfSegments',num_segments, 'Wavelength', wavelength, 'DetDist_mm', detector_distance, ...
        'PixelSize_mm', detector_pixelsize, 'FilenameValidMask',maskfilename, 'FilenameIntegMasks',integmaskfilename, ...
        nexus_args{:},'ConvertLocationArgs',convert_location_args);
        
    % Reintegrate silver behenate
    beamline.integrate_range_nexus(todo(1),'stopScanNr',todo(1),...
        'detector', detector,'outdirData',integdir,'basePath',datadir,'integMask',integmaskfilename);
end
%% Plot radial integration of Ag Behenate along with tabulated peak values

if Ifsaxs
    plotting.plot_radial_integ(sprintf('%s%s%d_%05d_00000_00000_integ.h5',integdir,userID,1,todo(1)),'PlotQ',1);
    numorders = numel(peaks2)-1;
    hold on
    plotting.vline(2*pi./(58.38./[1:numorders]));
    hold off
    legend({'data vs. 2\pin/(58.38 A^{-1})'})
else
    beamline.integrate_range_nexus(todo(1),'stopScanNr',todo(1),...
        'detector', detector,'outdirData',integdir,'basePath',datadir,'integMask',integmaskfilename);
    plotting.plot_radial_integ(sprintf('%s%s%d_%05d_00000_00000_integ.h5',integdir,userID,1,todo(1)),'PlotQ',1);
    hold on
    plotting.vline(2*pi./(58.38./[order_AgBe:11]));
    hold off
    legend({'data vs. 2\pin/(58.38 A^{-1})'})
end

%% step 5: radial integration & averaging of files -- 
%start here again if you merely want to integreat
%for fast measurements (i.e. scanning SAXS) start on several cn parallel
%adjust therefor integrate_range(scan_no_from,scan_no_to,scan_no_step) 
%and rund only step 0 and step 5
% save_format = '-v6';
%
% For a typical scanning SAXS experiment at 30 Hz it could be 1 computing
% node is enough to catch up
% save_format = '-v6';
starting_scan = min([todo scannr]);
fprintf('\n\nbeamline.integrate_range_nexus(%u,''detector'',''%s'',''outdirData'',''%s'',''basePath'',''%s'',''integMask'',''%s''',starting_scan, detector, integdir,datadir,integmaskfilename)
fprintf(');\n')

%%
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
