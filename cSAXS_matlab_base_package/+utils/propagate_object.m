% propagate_object.m

clear all   % Recommended if improfile will be used (there is a bug with cursor positioning otherwise)

addpath('..')
addpath('../ptycho/utils/')
import utils.*

%%% paths %%%
base_path='~/Data10/';
fig_path = '~/Data10/analysis/figures/';   % output directory for figures


%%% general settings %%%

% scan numbers
suffix = 'run_1_recons';
scan = [797];                 % scan number or range, ex.[464:794]

% profile
rangez = [-16e-3 16e-3];      % range (mm) for propagation profile
step_num = 200;               % number of steps for propagation profile
profile_type = 'amp';         % Either 'amp' or 'phase'.
% Propagate to one plane
prop_dis = -2.2e-3;           % chosen propagation distance for saving; leave empty for estimated distance

% Display option
disp = 'amp';               % Either 'hsv', 'amp', 'phase'
cmap = 'bone';              % colormap

% output
savedata=1;     % save output to file


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DO NOT EDIT THIS SECTION %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loading reconstruction

for kk=scan

    % get scan number
    scan =kk;
    
    % find filename 
    thisname=find_ptycho_filename(sprintf([base_path 'analysis/'], scan ), scan, [], suffix);
    
    % if file does not exist, wait for 10 seconds and try again
    while isempty(thisname)
        fprintf('Could not find any file. Waiting for 10 seconds...\n');
        pause(10)
        thisname=find_ptycho_filename(sprintf([base_path 'analysis/'], scan ), scan, [], suffix);
    end
    
    % if there are more than one file, just take the first one
    if iscell(thisname)
        thisname = thisname{1};
    end

    % get filename and load data
    [~, filename, ext]=fileparts(thisname);
    title_str = filename;
    filename = [filename ext];
    h = io.load_ptycho_recons(thisname);
    object = h.object;
    probe = h.probe;
    p = h.p;
    
    clear h;
    
    % get references from p
    lambda=p.lambda;            % wavelength [m]
    dis=p.z;                    % sample-detector distance [m]
    asize=p.asize(1);
    pixsize=p.dx_spec(1);       % pixel size of final reconstruction [m]
    
    % crop object borders
    ob_dims = size(object);
    %in square case
    %ob_dims = [min(size(object)) min(size(object))];
    %object = utils.crop_pad(object, ob_dims);

    apod=50;
    %mask2=fftshift(fract_hanning_pad(ob_dims(1),ob_dims(1)-asize+apod,ob_dims(1)-asize)); % square case
    mask2=fftshift(filt2d_pad(ob_dims,ob_dims-asize+apod,ob_dims-asize)); % 2D array
    ob=object.*mask2;

    scrsz = get(0,'ScreenSize');
    
    [ob_corr res]=remove_linearphase_v2(ob,ones(size(ob)),20);
    
    if isempty(prop_dis)
        prop_dis = utils.prop2focus(ob_corr, lambda, pixsize, 'd_start', 3.1e-3);
    end
    
    
    propdists = linspace(rangez(1),rangez(2),step_num);
    
    switch profile_type
        case 'amp'
            profile_func = @abs;
        case 'phase'
            profile_func = @angle;
    end
    
    for ii=1:step_num
        this_step=propdists(ii);
        prop_profile(:,ii) = sum(profile_func(prop_free_nf(ob, lambda, this_step, pixsize)),1);
    end
    
    figure(333);
    clf
    x = [-round(ob_dims(2)/2) round(ob_dims(2)/2)]*pixsize;
    % y = [-round(ob_dims(1)/2) round(ob_dims(1)/2)]*pixsize;
    imagesc(propdists*1e3,x*1e6,prop_profile)
    colormap(cmap)

    title_str = strrep(title_str, '_', ' ');
    title(['Axial propagation ' title_str]);
    xlabel('z [mm]')
    ylabel('x [\mum]')
    set(gcf,'Outerposition',[1 scrsz(4)-480 1000 480])

    % back-propagated image
    back = prop_free_nf(ob, lambda, prop_dis, pixsize);

    % Make plots
    back_sel=back(asize/2:size(back,1)-asize/2,asize/2:size(back,2)-asize/2);
    ob_sel=ob(asize/2:size(ob,1)-asize/2,asize/2:size(ob,2)-asize/2);
    ob_sel_dims=size(ob_sel);
    x_sel = [-round(ob_sel_dims(2)/2) round(ob_sel_dims(2)/2)]*pixsize;
    y_sel = [-round(ob_sel_dims(1)/2) round(ob_sel_dims(1)/2)]*pixsize;

    % propagated object
    figure(111);
    clf
    switch disp
        case 'hsv'
            imagesc(x_sel*1e6,y_sel*1e6,plotting.c2image(ob_sel));
        case 'amp'
            imagesc(x_sel*1e6,y_sel*1e6,abs(ob_sel));
            colormap(cmap)
        case 'phase'
            imagesc(x_sel*1e6,y_sel*1e6,angle(ob_sel));
            colormap(cmap)
    end
    
    axis xy equal tight;
    title(['object ' title_str]);
    xlabel('x [\mum]')
    ylabel('y [\mum]')
    colorbar
    set(gcf,'Outerposition',[501 1 500 500])


    % object without propagation
    figure(222);
    clf
    switch disp
        case 'hsv'
            imagesc(x_sel*1e6,y_sel*1e6,plotting.c2image(back_sel));
        case 'amp'
            imagesc(x_sel*1e6,y_sel*1e6,abs(back_sel));
            colormap(cmap)
        case 'phase'
            imagesc(x_sel*1e6,y_sel*1e6,angle(back_sel));
            colormap(cmap)
    end
    
    
    axis xy equal tight
    title(['object propagated to ' num2str(prop_dis*1e3) ' mm ' title_str]);
    xlabel('x [\mum]')
    ylabel('y [\mum]')
    colorbar
    set(gcf,'Outerposition',[1 1 500 500])

%% Save plots
if savedata 
    % define output directory
   [savedata_path, savedata_file] = fileparts(thisname);
   
   % define file name
   savedata_file_prop = [savedata_file sprintf('_obj_prop_%.2f_mm_test_', prop_dis*1e3)];
   savedata_file_prop = strrep(savedata_file_prop, '.', 'p');

   % save figures in fig_path
   if ~exist(fig_path,'dir')
       mkdir(fig_path)
   end
   
    %print('-f3','-djpeg', '-r300', fullfile(fig_path, [savedata_file '.jpg']))
    %print('-f2','-djpeg', '-r300', fullfile(fig_path, [savedata_file_prop '.jpg']))

    % save mat file 
    object=crop_pad(back, p.object_size);
    p.object_orig = p.object;
    p.object{1} = object;
    save(fullfile(savedata_path,[savedata_file_prop,'.mat']),'object','probe','p');

end


end


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
%   using the “cSAXS software package” developed by the CXS group,
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
