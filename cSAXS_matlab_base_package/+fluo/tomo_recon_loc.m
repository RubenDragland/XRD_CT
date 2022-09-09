%TOMO_RECON_LOC    this is a copy of tomo reconstruction code that was
%                  required to do tests of tomo fluo without creating a 
%                  separate branch in tomo package. Will be depricated 
%                  in later versions
%
% ** volume_in             aligned projections for reconstruction
% ** theta                 angles of projections
% ** alignment             par structure from aligned ptychography
%
% returns:
% ++ volume_out            reconstructed volume
%

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
function volume_out = tomo_recon_loc(volume_in, theta, alignment)
if rem(size(volume_in, 3), 2)
    disp('GPU tomo recon only works with an even number of projections. Removing the last projection.')
    volume_in = volume_in(:, :, 1:end-1);
end
par = alignment.par;
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Edit this section %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
par.usecircle = true;               % Use circle to mask out the corners of the tomogram
par.filter_type = 'ram-lak';        % FBP filter (ram-lak, hamming, ....)
par.freq_scale = 1;                 % Frequency cutoff
par.phase_unwrapping = 'slow_sino';

apodize = 0;                        % axial apodization 
radial_smooth_apodize = alignment.Npix(1)/100; % smooths the circle function (circulo), which delimits the reconstruction, this to avoid sharp edges in the  reconstruction

circulo = [];  % apodization function 
tomo_size = size(volume_in);
sinogram = zeros(tomo_size); %zeros([tomo_size(2) tomo_size(3) tomo_size(1)]);

N_slice = size(volume_in,1);
for slice_num = 1:N_slice
    clear sino
    sino = squeeze(volume_in(slice_num,:,:));
    sinogram(slice_num,:,:) = sino;
end

[Ny_sino,Nx_sino,~] = size(sinogram);
CoR = [Ny_sino,Nx_sino]/2+[0,0.5]*strcmpi(par.phase_unwrapping,'none');
[cfg, vectors] = astra.ASTRA_initialize([alignment.Npix, alignment.Npix, Ny_sino],[Ny_sino,Nx_sino],theta,par.lamino_angle,par.tilt_angle,1,CoR); 

split = astra.ASTRA_find_optimal_split(cfg, length(par.GPU_list), 2);

if par.usecircle
    [~,circulo] = utils.apply_3D_apodization(ones(alignment.Npix), apodize, 0, radial_smooth_apodize); 
end

tomograms = {};
for ii = 1:2
    tomograms{ii} = tomo.FBP_zsplit(sinogram, cfg, vectors, split,'valid_angles', [],...
        'GPU', par.GPU_list,'filter',par.filter_type, 'filter_value',par.freq_scale,...
        'use_derivative', strcmpi(par.phase_unwrapping, 'none'), 'mask', circulo);
end


volume_out = ((tomograms{1} + tomograms{2})/2);
volume_out = utils.crop_pad(volume_out, [97 97]);

end