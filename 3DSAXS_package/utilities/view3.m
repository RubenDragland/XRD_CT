% function view3(I,coloraxis,viewpar)
% Function to visualize 3D data slices interactively
% Inputs
% I                 3D array
% coloraxis         Axis for the colorbar, e.g. = [-1 1] or 'auto'
% viewpar           Structure with other parameters for visualization
% viewpar.which     Vector with position where cuts are displayed, default [round(sizey/2) round(sx/2) round(sz/2)]
% viewpar.interactive   If (= true, default) you can click the slices and
%                       see the corresponding perpendicular slices corresponding to the point you
%                       clicked. If (= false) then it just plots once.

%*------------------------------------------------------------------------*
%|                                                                        |
%|  Except where otherwise noted, this work is licensed under a           |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0             |
%|  International (CC BY-NC-SA 4.0) license.                              |
%|                                                                        |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)     |
%|                                                                        |
%|      Author: CXS group, PSI                                            |
%*------------------------------------------------------------------------*
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function view3(I,coloraxis,viewpar)

colors = [1 0 0; 0 0.7 0; 1 0 1];

[sy sx sz] = size(I);
if ~exist('viewpar')
    viewpar.nothing = [];
end
if ~isfield(viewpar,'which')
    which = [round(sy/2) round(sx/2) round(sz/2)];
else
    which = viewpar.which;
end

if ~isfield(viewpar,'interactive')
    interactive = true;
else
    interactive = viewpar.interactive;
end

reply = 1;
while (reply==1)
    if ~interactive
        reply = 0;
    end
    colormap jet
    figure(1);
    subplot(2,2,1)
    imagesc(I(:,:,which(3)));
    axis xy equal tight
    colorbar
    title(['XY [1], Z = ' num2str(which(3))])
    set(gca,'XColor',colors(3,:),'YColor',colors(3,:))
    xlabel('X')
    ylabel('Y')
    hold on,
    plot([1 sx],[1 1]*which(1),'--','Color',colors(1,:))
    plot([1 1]*which(2),[1 sy],'--','Color',colors(2,:))
    hold off,
    h1 = gca;
    
    subplot(2,2,2)
    imagesc(squeeze(I(:,which(2),:)));
    axis xy equal tight
    colorbar
    title(['YZ [2], X = ' num2str(which(2))])
    set(gca,'XColor',colors(2,:),'YColor',colors(2,:))
    xlabel('Z')
    ylabel('Y')
    hold on,
    plot([1 sx],[1 1]*which(1),'--','Color',colors(1,:))
    plot([1 1]*which(3),[1 sy],'--','Color',colors(3,:))
    hold off,
    h2 = gca;
    
    subplot(2,2,3)
    imagesc(flipud(rot90(squeeze(I(which(1),:,:)),1)));
    axis xy equal tight
    colorbar
    title(['XZ [3], Y = ' num2str(which(1))])
    xlabel('X')
    ylabel('Z')
    set(gca,'XColor',colors(1,:),'YColor',colors(1,:))
    hold on,
    plot([1 sx],[1 1]*which(3),'--','Color',colors(3,:))
    plot([1 1]*which(2),[1 sz],'--','Color',colors(2,:))
    hold off,
    h3 = gca;
    
    if exist('coloraxis','var')
        subplot(2,2,1)
        caxis(coloraxis)
        subplot(2,2,2)
        caxis(coloraxis)
        subplot(2,2,3)
        caxis(coloraxis)
    end
    %reply = input('Which plot will you pick? [1 2 3]   ');
    if reply
        selected = ginput(1),
    else
        selected = [];
    end
    if ~isempty(selected)
        hsel = gca;
        switch hsel
            case h1
                which(2) = round(selected(1));
                which(1) = round(selected(2));
            case h2
                which(3) = round(selected(1));
                which(1) = round(selected(2));
            case h3
                which(2) = round(selected(1));
                which(3) = round(selected(2));
            otherwise
                break;
        end
    else
        reply = 0;
    end
        
end

