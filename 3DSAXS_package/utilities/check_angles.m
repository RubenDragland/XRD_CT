% function [output_tomo] = check_angles(filename,viewpar)
% Function to visualize output of optimization using SH,
% visualize 3D data slices interactively showing 3D orientation in color,mapped on one hemisphere
% as well as in the selected voxel the plot of 3D reciprocal space map
% modelled with spherical harmonics
% Inputs
% filename              path and name of output file from SH optimization
%                       containing s and p.lines and p.points_per line
% viewpar               Structure with other parameters for visualization
% viewpar.which         Vector with position where cuts are displayed, default [round(sizey/2) round(sx/2) round(sz/2)]
% viewpar.interactive   If (= true, default) you can click the slices and
%                       see the corresponding perpendicular slices corresponding to the point you
%                       clicked. If (= false) then it just plots once.
% viewpar.mask3D        3D binary array, puts value of color to black if
%                       mask value = 0
% viewpar.symnint       If (=true) 3 puts value of color to the symmetric intensity (first coefficient of a_out1)
%                       default is false
% viewpar.dego          If (=true) puts value of color to degree of
%                       orientation calcuated from aout_1
%                       default is false
% viewpar.plotradius    If false (=default) the spherical harmonic function is shown on a sphere, 
%                       if true the radius is proportional to the function
%
% Output                structure containing the reshaped a_tomo, theta_tomo and phi_tomo
% uses the functions coloring_3D_plot, plot_spherical_harmonic_squared
% Function to visualize 3D data slices
%*-------------------------------------------------------------------------------------*
%|                                                                                                           |
%|  Except where otherwise noted, this work is licensed under a            |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0        |
%|  International (CC BY-NC-SA 4.0) license.                                         |
%|                                                                                                           |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)       |
%|                                                                                                           |
%|      Author: CXS group, PSI                                                                |
%*------------------------------------------------------------------------------------*
% Version 5.0
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
function [output_tomo] = check_angles(filename,viewpar)


load(filename)

% compatibility fix
theta_out1 = s.theta.data;
phi_out1 = s.phi.data;
a_out1 = [];
for ii = 1:numel(s.a)
    a_out1 = [a_out1 s.a(ii).data];
end

output_tomo.a_tomo = reshape(a_out1,p.points_per_line,p.lines,p.lines,[]);
output_tomo.theta_tomo = reshape(theta_out1,p.points_per_line,p.lines,p.lines);
output_tomo.phi_tomo = reshape(phi_out1,p.points_per_line,p.lines,p.lines);
output_tomo.p = p;
a_tomo = output_tomo.a_tomo;
phi_tomo = output_tomo.phi_tomo;
theta_tomo = output_tomo.theta_tomo;
colors = [1 0 0; 0 0.7 0; 1 0 1];

[sy sx sz] = size(theta_tomo);

if size(a_tomo,4) == 4;
    l = [0,2,4,6];
    m = [0,0,0,0];
elseif size(a_tomo,4) == 1;
    l = [0];
    m = [0];
elseif size(a_tomo,4) == 2;
    l = [0,2];
    m = [0,0];
elseif size(a_tomo,4) == 3;
    l = [0,2,4];
    m = [0,0,0];
end

if exist('viewpar') == 0;
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

if ~isfield(viewpar,'plotradius') 
    plotradius = 0;    
else
    plotradius = viewpar.plotradius;
end

if isfield(viewpar,'symint') && ~isempty(viewpar.symint) && ~viewpar.symint==0
    mask_3D = a_tomo(:,:,:,1)./max(max(max(a_tomo(:,:,:,1))));
end

if isfield(viewpar,'dego') && ~isempty(viewpar.dego) && ~viewpar.dego==0
    if isfield(viewpar,'symint')
        fprintf('you have chosen viewpar.symint and viewpar.dego only the later will be applied \n');
    end
    if size(a_tomo,4) == 4
        degorientation = (abs(a_tomo(:,:,:,2).^2) + abs(a_tomo(:,:,:,3).^2) + abs(a_tomo(:,:,:,4).^2))./(abs(a_tomo(:,:,:,1).^2) + abs(a_tomo(:,:,:,2).^2) + abs(a_tomo(:,:,:,3).^2) + abs(a_tomo(:,:,:,4).^2));
    elseif size(a_tomo,4) == 3
        degorientation = (abs(a_tomo(:,:,:,2).^2) + abs(a_tomo(:,:,:,3).^2))./(abs(a_tomo(:,:,:,1).^2) + abs(a_tomo(:,:,:,2).^2) + abs(a_tomo(:,:,:,3).^2));
    else
        fprintf('adjust script for right order \n');
    end
    mask_3D = degorientation./max(degorientation(:));
end

if exist('mask_3D') == 0;
    mask_3D = ones(size(theta_tomo));
end

if isfield(viewpar,'mask3D') && ~isempty(viewpar.mask3D)
    mask_3D = mask_3D .*viewpar.mask3D;
end

reply = 1;
while (reply==1)
    if ~interactive
        reply = 0;
    end
    
    figure(2);
    
    subplot(2,2,1)
    
    coloring_3D_plot(theta_tomo(:,:,which(3)),phi_tomo(:,:,which(3)),mask_3D(:,:,which(3)),0,1);
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
    coloring_3D_plot(squeeze(theta_tomo(:,which(2),:)),squeeze(phi_tomo(:,which(2),:)),squeeze(mask_3D(:,which(2),:)),0,1);
    
    
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
    coloring_3D_plot(squeeze(theta_tomo(which(1),:,:)),squeeze(phi_tomo(which(1),:,:)),squeeze(mask_3D(which(1),:,:)),0,1);
    title(['XZ [3], Y = ' num2str(which(1))])
    xlabel('X')
    ylabel('Z')
    set(gca,'XColor',colors(1,:),'YColor',colors(1,:))
    hold on,
    plot([1 sx],[1 1]*which(3),'--','Color',colors(3,:))
    plot([1 1]*which(2),[1 sz],'--','Color',colors(2,:))
    hold off,
    h3 = gca;
    
    subplot(2,2,4)
    %title([output.p.add_name])
    plot_spherical_harmonic_squared([theta_tomo(which(1),which(2),which(3)),phi_tomo(which(1),which(2),which(3))],l,m,[squeeze(a_tomo(which(1),which(2),which(3),:))],plotradius);
    colormap jet
    
    
    %reply = input('Which plot will you pick? [1 2 3]   ');
    selected = ginput(1),
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

