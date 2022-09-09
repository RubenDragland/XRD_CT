%% [tomotensor, B_segs, error_overall] = tensor_reconstruct(projection,  num_iter, if_show, initial_guess,subset, parallel_scattering)
% Reconstructs the tensor model from SAXS projections.
%
% Inputs:
%   projection  structure containing data of SAXS projections.
%   num_iter    number of iteratoins, usually ~10000 or less.
%
% Above is required for user input.
%   if_show     display the model during reconstruction, =0 by default.
%   initial_guess  starting model for the reconstruction, can be used to
%   continue from previous result, empty by default.
%   subset      subset of projection angles to be used for reconstruction,
%   =[] for all angles (default).
%   parallel_scattering     =1 to assume scattering parallel to the
%   structural model (required for combination with SASTT). =0 by default.

% Copyright 2017 Zirui Gao
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [tomotensor, B_segs, error_overall] = tensor_reconstruct(projection,  num_iter, if_show, initial_guess,subset, parallel_scattering)
% reconstruct a tensor model from all projections
if nargin<6
    parallel_scattering = 0;
end
% =0 to set at default perpendicular scattering.

% Calculate the projection B matrix
% The details about this step can be found in 'calc_B_segs.m'.
% this is like the magnetization vector from NMR
[B_segs] = calc_B_segs(projection,parallel_scattering);

% this is the regularization factor
ratio = 0.01; %(1/numprojs);   % ratio of error correction, approximately 0.01, optimal value can be selected from testing with error list

if nargin < 5
    subset = [];
end




% Modify until here.
numprojs = length(projection);
err=zeros(num_iter);
which_proj=zeros(num_iter);
samp = projection(1).data;
if nargin<3
    if_show=0;
end

if (nargin<4)||isempty(initial_guess)
    tomotensor = zeros(size(samp,1), size(samp,2), size(samp,2), 6);
else
    tomotensor = initial_guess;
end

if isempty(subset)
    subset = [1:numprojs];
end
proj_shifted=cell(length(projection),1);
T_target=cell(length(projection),1);


fprintf('Shifting projections for alignment\n');
tic;
for ii = 1:length(projection)
    samp=projection(ii).data;
    projection(ii).Rot_exp = double(projection(ii).Rot_exp);
    for kk = 1:size(projection(ii).data, 3)
        if isfield(projection(1),'window_mask')
            proj_shifted{ii}(kk,:,:) = double(utils.shiftwrapbilinear(squeeze(projection(ii).data(:,:,kk)).*projection(ii).window_mask, ...
                -projection(ii).dy, -projection(ii).dx));
        else
            proj_shifted{ii}(kk,:,:) = double(utils.shiftwrapbilinear(squeeze(projection(ii).data(:,:,kk)), ...
                -projection(ii).dy, -projection(ii).dx));
        end        
    end
    if isfield(projection(1),'window_mask')
        projection(ii).window_mask=double(imtranslate(projection(ii).window_mask, ...
            round([projection(ii).dx, projection(ii).dy])));
    end
    T_target{ii}=zeros(6,size(samp,1),size(samp,2));
    for x=1:size(samp,1)
        for y=1:size(samp,2)
            T_target{ii}(:,x,y)=squeeze(B_segs(ii,:,:))'*squeeze(proj_shifted{ii}(:,x,y));
        end
    end
            
end
toc;
%%
tic;
rng shuffle;
for iter = 1:num_iter
    % choses a random projection
    i = randi(size(subset,2));
    i = subset(i);
    
    B = squeeze(B_segs(i,:,:));
    
    % find the random projection
    prawt = T_target{i};
    % correct the segments
    %prawt = prawt.*projection(i).window_mask_shifted;
    
    % pemute as the input for the function = [segments Y X]
    % gets the dimensions [Y, X]
    projsize = size(prawt);
    psimt = proj_simul_tensor(tomotensor, projection(i).Rot_exp, B'*B, projsize);
    
    diff = (prawt-psimt);
    err(iter)=rms(diff(:));
    which_proj(iter)=i;
    
    if isfield(projection(1),'window_mask')
        for s=1:6
            diff(s,:,:)=squeeze(diff(s,:,:)).*projection(i).window_mask;
        end
    end
    
    
    diff = diff.*ratio;
    
    % tomotensor has dimensions [Y X Z seg]
    tomotensor = proj_tensor_correct(tomotensor,projection(i).Rot_exp, diff, projsize);
    if (mod(iter,100)==0)
        fprintf('Iteration %d \n', iter);
    end
    if (mod(iter,100)==0)&&(if_show==1)
        proj_rawt = proj_shifted{i};
        projsize = size(proj_rawt);
        proj_simt = proj_simul_tensor(tomotensor, projection(i).Rot_exp, B, projsize);
        if isfield(projection(1),'window_mask')
            window_mask = zeros(projsize);
            for ii = 1:size(window_mask, 1)
                window_mask(ii, :,:) = projection(i).window_mask;
            end
        else
            window_mask = ones(projsize);
        end
        figure(56)
        subplot(2,3,1)
%         imagesc(squeeze(mean(prawt,1)));
        fourier_figure_fast(proj_rawt.*window_mask,projection(1),gcf);
        ax=findall(gcf,'type','axes');
        axes(ax(2));
        title(sprintf('REF: rotX%.0f, rotY%.0f', projection(i).rot_x, projection(i).rot_y));
        subplot(2,3,2)
        fourier_figure_fast(proj_simt.*window_mask,projection(1),gcf);
        ax=findall(gcf,'type','axes');
        axes(ax(2));
        title(sprintf('OPT: rotX%.0f, rotY%.0f', projection(i).rot_x, projection(i).rot_y));

        
        diff_plot = (proj_rawt-proj_simt).*window_mask;
        subplot(2,3,3)
        imagesc(squeeze(mean(diff_plot, 1)));
        colorbar
        title(sprintf('It. %d: difference', iter));
        axis xy tight equal
                
        subplot(2,3,[4:6])
%         errorbar(iter, mean(err(iter-100+1:iter)), std(err(iter-100+1:iter)), 'ok--');    
        plot(err(1:iter));
%         hold on
        xlabel('Iteration')
        ylabel('Mean difference')
%         set(gca, 'yscale', 'log' )
%         title(sprintf('RotX = %.02f, Roty = %.02f', rotx, roty));
%         axis xy tight equal
        grid on
        drawnow
        
    end
    
end
toc;

% Calculate overall error
if nargout>2
    err_all=zeros(size(subset,2),1);
    for ii=1:size(subset,2)
        i=subset(ii);
        proj_rawt = proj_shifted{i};
        projsize = size(proj_rawt);
        B = squeeze(B_segs(i,:,:));
        proj_simt = proj_simul_tensor(tomotensor, projection(i).Rot_exp, B, projsize);
        diff_plot = (proj_rawt-proj_simt);
        err_all(ii)=rms(diff_plot(:));
    end
    error_overall=rms(err_all(:));
    fprintf('Overall error is %f\n',error_overall);
end

end