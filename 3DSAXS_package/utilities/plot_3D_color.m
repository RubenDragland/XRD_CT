% fig_save = plot_3D_color(s, slice_nr)
% Plot 3D orientation of tensor tomogram
%
% Inputs:
% ** s          Tensor structure
% ** slice_nr   Number of slice to display
%
% Output:
% ++ fig_save   Figure handle

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

function fig_save = plot_3D_color(s, slice_nr)

mask3D = s.mask3D;
s.phi.data = mod(s.phi.data, 2*pi);
s.theta.data = mod(s.theta.data, 2*pi);
[ny nx nz] = size(s.phi.data);

if isempty(slice_nr)
    slice_nr_x = round(nx/2);
    slice_nr_y = round(ny/2);
    slice_nr_z = round(nz/2);
else
    slice_nr_x = slice_nr;
    slice_nr_y = slice_nr;
    slice_nr_z = slice_nr;
end

% plot the results
fig_save = figure(45);
hold on

subplot(2, 2, 1)
coloring_3D_plot(squeeze(s.theta.data(:,:,slice_nr_z)), ...
    squeeze(s.phi.data(:,:,slice_nr_z)), squeeze(mask3D(:,:,slice_nr_z)), 0, 1);
title(sprintf('Slice %d', slice_nr_z))
xlabel('x')
ylabel('y')

subplot(2, 2, 2)
coloring_3D_plot(squeeze(s.theta.data(:,slice_nr_x,:)), ...
    squeeze(s.phi.data(:,slice_nr_x,:)), squeeze(mask3D(:,slice_nr_x,:)),0, 1);
title(sprintf('Slice %d', slice_nr_x))
xlabel('z')
ylabel('y')

subplot(2, 2, 3)
coloring_3D_plot((squeeze(s.theta.data(slice_nr_y,:,:)))', ...
    (squeeze(s.phi.data(slice_nr_y,:,:)))', squeeze(mask3D(slice_nr_y,:,:))',0, 1);
title(sprintf('Slice %d', slice_nr_y))
xlabel('x')
ylabel('z')

subplot(2, 2, 4)
coloring_3D_plot(0, 0, 0, 1, 0);
axis off
drawnow

