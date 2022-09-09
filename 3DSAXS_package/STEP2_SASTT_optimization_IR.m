% WIP: Optimization with the Ir method

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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent = cd ;
base_path = '/Data sets/' ; %'~/Data10'; % = ~/Data10 for online analysis, provide the path for offline analysis Ex: '/das/work/p16/p16649/'
base_path = [parent base_path] ; % SAFE WITH FULL PATH
sample_name = 'SASTT_carbon_knot_aligned_ASTRA_corrected' ; % name given in the saxs_caller_template
add_name = ''; % 'online'; i.e.       % additional name the optimizations: = [ ] if not needed
use_ASTRA = 1; % =1 if you already have results from ASTRA alignment,
               % =0 and will do tomo reconstruction and alignment from the
               % beginning.
                 
if_show = 1;   % =1 to show the steps in reconstruction.
               % =0 to not show any figure.
               
                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
% filename = fullfile(base_path,sprintf('analysis/SASTT/%s/SASTT_%s_aligned_ASTRA.mat', ...
%     sample_name, sample_name));
%     load(filename);

%% Tomo reconstruction and Projection Registration (or load directly from ASTRA)
if use_ASTRA==0
    %filename = fullfile(base_path,sprintf('analysis/SASTT/%s/projection_data/SASTT_%s.mat', ...
    %    sample_name, sample_name));
    filename = sprintf('%s%s.mat', base_path, sample_name);
    fprintf('Reconstructing tomogram by slices\n');
    [tomo,projection]=tomo_reconstruct(filename,if_show);
    fprintf('Registering projections to tomogram\n');
    projection_aligned=register_projections(projection,tomo);
else
    %filename = fullfile(base_path,sprintf('analysis/SASTT/%s/projection_data/SASTT_%s_aligned_ASTRA.mat', ...
    %sample_name, sample_name));
    filename = sprintf('%s%s.mat', base_path, sample_name);
    load(filename);
    projection_aligned=projection;
end

params_IRTT.results = fullfile(base_path, sprintf('analysis/SASTT/%s/IRTT/%s/results/', sample_name, add_name));
if ~exist(params_IRTT.results, 'file')
    mkdir(params_IRTT.results);
end
params_IRTT.figures = fullfile(base_path, sprintf('analysis/SASTT/%s/IRTT/%s/figures/', sample_name, add_name));
if ~exist(params_IRTT.figures, 'file')
    mkdir(params_IRTT.figures);
end

%% Tensor reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct tensor model from projection of all segments.
% Set a lower iteration first to check the result, if it's correct then
% run more iterations to minimize the error.

num_iter = 10000;
% Typical reconstruction would need 1000 iterations to see the rough
% model, and 10000 iterations for final result.

continue_recon = 0; 
%   =1 to continue the iterations from previous result.
%!! =0 will OVERWRITE previous reconstruction results.
            

subset= [];
% which subset of projections to use in the reconstruction, 
% =[] to use all projections.

parallel_scattering = 0; 
% If to reconstruct scattering parallel to perpendicular to the tensor.
% =0 by default.
% =1 is required if used for initial guess of SH.
                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tensor_filename = fullfile(sprintf('%sSASTT_%s_tomotensor_%s.mat', ...
    params_IRTT.results,sample_name, add_name));

if ~continue_recon
    [tomotensor, B_segs,error_overall] = tensor_reconstruct(projection_aligned, ...
        num_iter, if_show, [], subset,parallel_scattering);
else
    load(params_IRTT.results);
    [tomotensor, B_segs,error_overall] = tensor_reconstruct(projection_aligned, ...
        num_iter, if_show, tomotensor, subset,parallel_scattering);
end

filename = sprintf('%sIRTT_%s_q_%d-%d_%s', params_IRTT.figures, sample_name, projection(1).par.r_sum{1}(1),projection(1).par.r_sum{1}(end), add_name);
        print(gcf, filename, '-dpng','-r300');
% The format of tomotensor output is a  [Y, X, Z, T] array of a symmetric
% tensor for each voxel in the model.
% For the last dimension T the elements are arranged as (Txx,Tyy,Tzz,Txy,Txz,Tyz)

%% Save tensor model (will overwrite)

save(tensor_filename,'tomotensor','B_segs','error_overall','parallel_scattering');
fprintf("Saved %s\n",tensor_filename);

%% Tensor model output (eigen vector map)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visulizing the tensor directly is not easy, therefore here it solves a
% certain eigen vector of the tensor and visualize the vector map.

% This is done by writing into vtk files and load with Paraview.

which_eig = 1;
% Which eigenvector to solve out. =1 normally.

% The threshold for the eigenvalue to remove background.
threshold_low = 0;   
threshold_high = inf; 

filename_add = 'thres1'; % To add to the output filename, e.g. threshold.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the eigen vectors
eig_vecs = solve_eig(tomotensor, which_eig, threshold_low,threshold_high);

% Write the vtk files

% This part is only for visualizing in Paraview. Don't need to change.
resol=projection(1).par.x_scale*1000;
vtk_filename = fullfile(sprintf('%sSASTT_%s_Eig%d_%s.vtk', ...
    params_IRTT.results, sample_name, which_eig, filename_add));
vtkwrite(vtk_filename,resol,'VECTORS',sample_name,eig_vecs);

fprintf("Saved %s\n",vtk_filename);
