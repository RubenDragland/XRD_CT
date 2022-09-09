%%COMPARE_COLORMAPS plots your data or sample data with various colormaps
%
% *optional*
% inArr             2D input array
%
% EXAMPLES:
%       % dummy data
%       plotting.colormaps.compare_colormaps();
%       
%       % 2D user data 
%       plotting.colormaps.compare_colormaps(rand(256, 256));
%
% see also: plotting.colormaps.create_colormap

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
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

function compare_colormaps(varargin)

if nargin>0
    inArr = varargin{1};
    sz = size(inArr);
    if numel(sz)>3
        error('Please provide a 2D array.');
    end
else 
    inArr = [];
end

close all;

%% ignore the following files:
cmapIgnore = ['compare_colormaps.m', 'create_colormap.m'];

%% some Matlab colormaps
cmapsMatlab = {'bone', 'parula', 'jet'};

cPath = strsplit(mfilename('fullpath'), '/');
cPath = strjoin(cPath(1:end-1),'/');

%% get the cmaps
cmapList = [];
res = dir(fullfile(cPath, '*.m'));
for ii=1:length(res)
    if ~contains(cmapIgnore, res(ii).name)
        cmapList{end+1} = ['plotting.colormaps.' res(ii).name(1:end-2)];
    end
end
% add additional MATLAB cmaps
for ii=1:length(cmapsMatlab)
    cmapList{end+1} = cmapsMatlab{ii};
end

%% load data or plot user data
if ~isempty(inArr)
    img = repmat(inArr(:,:,1), [1 1 length(cmapList)]);
    figure();
    plotting.imagesc3D(img);
    ax_controller = gca;
    propListener = addlistener(ax_controller.slider_handle,'Value','PostSet',@(src,evnt)change_cmap(ax_controller, cmapList));
    set(ax_controller.slider_handle,'Value', 2);
    set(ax_controller.slider_handle,'Value', 1);
    
else
    
    controllerID = 5;
    
    
    %% load the sample data
    img{1} = imread(fullfile(cPath, 'private', 'ML512.jpg'));
    img{1} = img{1}(:,:,1);
    
    img{2} = load(fullfile(cPath, 'private', 'probe_PSI.mat'));
    img{2} = log10(abs(fftshift(fft2(img{2}.probe))));
    
    img{3} = imread(fullfile(cPath, 'private', 'chip_phantom.png'));
    img{3} = double(img{3}(:,:,1));
    
    img{4} = io.image_read(fullfile(cPath, 'private', 'ptycho_00000.h5'), 'H5Location', '/eh5/images');
    img{4} = log10(img{4}.data(420:820,50:450));
    
    img{5} = io.image_read(fullfile(cPath, 'private', 'SAXS_00051.cbf'));
    img{5} = log10(img{5}.data(400:1000,450:1050));
    
    [g1, g2] = utils.get_grid(256, 1);
    img{6} = angle(fftshift(exp(1j.*(g1.^2 + g2.^2)./256)));
    
    
    %% plot
    img{controllerID} = repmat(img{controllerID}, [1 1 length(cmapList)]);
    ax = {};
    figure('units','normalized','outerposition',[0 0.1 1 0.9]);
    for ii=1:6
        subplot(2,3,ii);
        if ii==controllerID
            ax_controller = gca;
            plotting.imagesc3D(img{ii}); axis equal tight;
            colorbar;
        else
            imagesc(img{ii}(:,:,1)); axis equal tight;
            colorbar;
            ax{end+1} = gca;
            ax{end}.addprop('img');
            ax{end}.img = img{ii};
        end
        
    end
    %% add listener
    propListener = addlistener(ax_controller.slider_handle,'Value','PostSet',@(src,evnt)change_cmap(ax_controller, cmapList, ax));
    set(ax_controller.slider_handle,'Value', 2);
    set(ax_controller.slider_handle,'Value', 1);
    
    
end
end

function change_cmap(ax, cmapList, varargin)
if ~isempty(varargin)
    ax1 = varargin{1};
    multiplot = true;
else
    multiplot = false;
end

try
    cmapFunc = str2func(cmapList{max(min(length(cmapList),round(ax.slider_handle.Value)),1)});
    if multiplot
        for ii=1:length(ax1)
            ax1{ii}.Colormap = cmapFunc();
        end
    end
    ax.Colormap = cmapFunc();
    plotting.suptitle(cmapList{max(min(length(cmapList),round(ax.slider_handle.Value)),1)}, 'Interpreter', 'none');
    
catch
end
end


