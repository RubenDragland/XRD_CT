% Calculate the degree of orieantation and export to paraview


%% Step 4.1 load the data and check the optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parent = cd ;
base_path = '/Data sets/' ; 
base_path = [parent base_path] ; % RSD: SAFE WITH FULL PATH
sample_name = 'SASTT_carbon_knot_aligned_ASTRA_corrected';
add_name = '';%'ID';       % additional name the optimizations: = [ ] if not needed
which_coeff = 1;                % which coefficient to show a0 = 1, a2 = 2, a4 = 3, a6 = 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath utilities/

%filename = fullfile(base_path,sprintf('analysis/SASTT/%s/projection_data/SASTT_%s_aligned_ASTRA.mat', sample_name, sample_name));
filename = sprintf('%s%s.mat', base_path, sample_name); %RSD: CURRENTLY WHAT IS NEEDED.
orig_data = load(filename);
% load the optimization results from SH
filename = fullfile(base_path, sprintf('analysis/SASTT/%s/SH/%s/results/', sample_name, add_name));
        filename = fullfile(filename, sprintf('result_%s_q%d-%d_all_again_%s.mat', sample_name, orig_data.projection(1).par.r_sum{1}(1), orig_data.projection(1).par.r_sum{1}(end), add_name));
clear orig_data
load(filename);

reshaped_coefficients = reshape(s.a(which_coeff).data, size(s.mask3D) ) ;
% to visualize the coefficients
%view3(s.a(which_coeff).data.*s.mask3D);  % RSD: ORIGINAL SCRIPT HAS FORGOTTEN TO RESHAPE
view3(reshaped_coefficients.*s.mask3D);

%% STEP 4.2 calculate degree of orientation

degree_orientation = zeros(size(s.theta.data)); 
%degree_orientation = zeros(size(s.mask3D));
num = degree_orientation;
den = degree_orientation;

for ii = 1:numel(s.a)
    if (s.a(ii).l ~= 0)||(s.a(ii).m ~= 0)
        num = num + abs(s.a(ii).data).^2;
    end
    den = den + abs(s.a(ii).data).^2;
end

degree_orientation = reshape( (num./den) , size( s.mask3D) ) ; %RSD: s.theta.data is 1D while we wish to display 3D
%degree_orientation = (num./den).*s.mask3D;
degree_orientation = degree_orientation.*s.mask3D;

view3(degree_orientation.*s.mask3D); % RSD: WHY IS THE MASK APPLIED TWICE?

%% Step 4.3 check the result of the degree of orientation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
apply_threshold = 0; % = 0 for applying no threeshowld to the degree of orientation
threshold_max = 0.5; % max degree of orientation
threshold_min = 0; % min degree of orientation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%apply thresholds to the degree of orientation
if (apply_threshold)
    degorientation_threshold = degree_orientation;
    degorientation_threshold(degorientation_threshold > threshold_max) = 0;
    degorientation_threshold(degorientation_threshold < threshold_min) = 0;
    view3(degorientation_threshold.*s.mask3D);
    degree_orientation = degorientation_threshold;
    clear degorientation_threshold degorientation_threshold_tomo
end


%% Step 4.4 Export data for paraview
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This will export vtk files of orientation vectors to be visualized in
% paraview.
magnitude_symIntensity = 0; % If to use symmetric intensity a0 as vector magnitude.
% This should be =1 in most cases.
magnitude_degOrientation = 1; % If to scale the vectors with degree of orientation.
% This depends on the dataset.
thres_low = 2.5;
thres_high = inf;
% Thresholds for magnitude of final output vectors, can be checked in paraview.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Strange path to vtkwrite should be changed   MGS
addpath utilities/IRTT/utilities/

% Look for the index of symmetrical component
for ii = 1:numel(s.a)
    if (s.a(ii).l == 0)&&(s.a(ii).m == 0)
        ind_sym = ii;
        break
    end
end

% transform into vector data
xx_tomo = sin(s.theta.data).*cos(s.phi.data);
yy_tomo = sin(s.theta.data).*sin(s.phi.data);
zz_tomo = cos(s.theta.data);

if magnitude_symIntensity&&magnitude_degOrientation
    error('You have selected both outputs of symIntensity and degOrientations. Choose only one of them.')
end

if magnitude_symIntensity
    xx_tomo = xx_tomo.* s.a(ind_sym).data;
    yy_tomo = yy_tomo.* s.a(ind_sym).data;
    zz_tomo = zz_tomo.* s.a(ind_sym).data;
    extra_title = '_symIntensity';
end

if magnitude_degOrientation
    xx_tomo = xx_tomo.* degree_orientation;
    yy_tomo = yy_tomo.* degree_orientation;
    zz_tomo = zz_tomo.* degree_orientation;
    extra_title = '_deg_Orientation';
end

data_paraview = cat(4, yy_tomo, xx_tomo, zz_tomo);
data_paraview (isnan(data_paraview)) = 0;

%%% export as vtk file

% Apply thresholds to the norm
data_paraview_temp = data_paraview;
for ii=1:p.ny
    for jj=1:p.nx
        for kk=1:p.nz
            if (norm(squeeze(data_paraview(ii,jj,kk,:)))<thres_low)||(norm(squeeze(data_paraview(ii,jj,kk,:)))>thres_high)
                data_paraview(ii,jj,kk,:)=0;
            end
        end
    end
end

vtk_filename = sprintf('%sanalysis/SASTT/%s/SH/%s/optimization_output/SH_paraview_data_%s_%s%s.vtk', base_path, sample_name,...
    add_name, sample_name, add_name, extra_title);
vtkwrite(vtk_filename, 1, 'VECTORS',sample_name,data_paraview);
fprintf('Saved %s\n',vtk_filename);

%% alternative to export as raw file
%concatenate matrices
data_paraview = cat(4, xx_tomo, yy_tomo, zz_tomo);
data_paraview = permute(data_paraview, [4 2 1 3]);

% save data
filename = fullfile(p.optimization_output_path, sprintf ('paraview_data_%s_%s.raw', sample_name, add_name));
fprintf('Exporting data for ParaView in %s \n', filename);
fid = fopen(filename, 'w+');
fwrite(fid, data_paraview,'float');
fclose(fid);
fprintf('********************************\n');
fprintf('the tomogram size is: x = 0-%d\n', size(data_paraview, 2)-1)
fprintf('                      y = 0-%d\n', size(data_paraview, 3)-1)
fprintf('                      z = 0-%d\n', size(data_paraview, 4)-1)
fprintf('                scalars = %d\n', size(data_paraview, 1))

%%
%*-------------------------------------------------------------------------------------*
%|?? ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????               |
%|?? Except where otherwise noted, this work is licensed under a???? ??????????????              |
%|?? Creative Commons Attribution-NonCommercial-ShareAlike 4.0 ????????????                   |
%|?? International (CC BY-NC-SA 4.0) license. ????????????????????????????????????????????????????????              |
%| ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????              |
%|?? Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)??????               |
%|???????????????????????????????????? ????????????????????????????????????????????????????????????????????????????????????????????????????????              |
%|?????????? Author: CXS group, PSI                                                         |
%*------------------------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If this code, or subfunctions or parts of it, is used for research in a
% ????publication or if it is fully or partially rewritten for another
% ????computing language the authors and institution should be acknowledged
%???? in written form and additionally you should cite:
% ????????M. Liebi, M. Georgiadis, A. Menzel, P. Schneider, J. Kohlbrecher,
% ????????O. Bunk, and M. Guizar-Sicairos, ???Nanostructure surveys of
%?????? ??macroscopic specimens by small-angle scattering tensor tomography,???
%?????? ??Nature 527, 349-352 (2015).???? (doi:10.1038/nature16056)
%
% A publication that focuses on describing features, or parameters, that
% ??????are already existing in the code should be first discussed with the
% ??????authors.
% ????
% This code and subroutines are part of a continuous development, they
%?????? are provided ???as they are??? without guarantees or liability on part
%?????? of PSI or the authors. It is the user responsibility to ensure its
% ??????proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%