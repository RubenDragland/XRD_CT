
%*-------------------------------------------------------------------------------------*
%|                                                                                     |
%|  Except where otherwise noted, this work is licensed under a                        |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0                          |
%|  International (CC BY-NC-SA 4.0) license.                                           |
%|                                                                                     |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)                  |
%|                                                                                     |
%|      Author: CXS group, PSI                                                         |
%*------------------------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If this code, or subfunctions or parts of it, is used for research in a
%   publication or if it is fully or partially rewritten for another
%   computing language the authors and institution should be acknowledged
%   in written form and additionally you should cite:
%     M. Liebi, M. Georgiadis, A. Menzel, P. Schneider, J. Kohlbrecher,
%     O. Bunk, and M. Guizar-Sicairos, “Nanostructure surveys of
%     macroscopic specimens by small-angle scattering tensor tomography,”
%     Nature 527, 349-352 (2015).   (doi:10.1038/nature16056)
% and
%     M. Liebi, M. Georgiadis, J. Kohlbrecher, M. Holler, J. Raabe, I.
%     Usov, A. Menzel, P. Schneider, O. Bunk and M. Guizar-Sicairos,
%     "Small-angle X-ray scattering tensor tomography: model of the
%     three-dimensional reciprocal-space map, reconstruction algorithm
%     and angular sampling requirements," Acta Cryst. A74, 12-24 (2018).
%     (doi:10.1107/S205327331701614X)
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%
% This code and subroutines are part of a continuous development, they
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function  [errormetric, penalty_error, regularization_coeff]=find_regularization_mu(projection, mu, it_max, p, s)

errormetric=[];                 % Arrays to store values
penalty_error=[];               % Arrays to store values
regularization_coeff=[];        % Arrays to store values
%end
filename = sprintf('%s/regularization/', p.optimization_output_path);
if ~exist(filename, 'file')
    mkdir(filename);
end

p.itmax = it_max;              % maximum number of iteration

% parameters for optimization of angles with the 4 coefficients
%%% Vector with spherical harmonic parameters
p.l = [0  2 4 6];  % Polar order
p.m = [0  0 0 0];  % Azimuthal order

p.find_orientation = 1;    % Optimize over the main orientation of the structure
p.init_guess = 0;          % in case we use initial guesses

p.opt_coeff = [1, 1, 1, 1];

p.find_coefficients = any(p.opt_coeff);   % Optimize over coefficients
p.avoid_wrapping = 1;      % avoid wrapping over 2Pi of the angle

p.regularization = 0;      % Sieves regularization on coeff
p.regularization_angle = 1;%regularization of angle

%parameters for the sieves regularization (blurring the gradient)
kernel3D=window3(5,5,5,@hamming);
p.kernel=kernel3D./sum(kernel3D(:)); % for normalization (sum equals 1)

% save the initial solution to be used at the start of each interation
s_init = s;

ticnow = tic;
for kk = 1:length(mu)
    p.regularization_angle_coeff = mu(kk);
    
    %optimize
    [p, s] = optimization.optimize_SH(projection, p, s_init);
    
    if p.regularization_angle_coeff == mu(end)
        filename = sprintf('%s/regularization/optimization_lcurve_%s_reg_%.0e', p.optimization_output_path, p.sample_name,...
            p.regularization_angle_coeff);
        print(gcf, filename,'-dpng','-r300');
        fprintf('Saving figure to: %s \n', filename)
    end
    
    % save results
    e = optimization.errorplot([]);
    filename = sprintf('%s/regularization/result_%s_q%d-%d_lcurve_%s_reg_%.0e.mat', p.optimization_output_path, p.sample_name, projection(1).par.r_sum{1}(1),  projection(1).par.r_sum{1}(end),...
        p.add_name, p.regularization_angle_coeff);
    save(filename,'s','e', 'p', '-v6');
    fprintf('Saving data to: %s \n', filename)
    
    %plot_3D_color(p, slice_nr, filename)
    fig_save = plot_3D_color(s, []);
    filename = sprintf('%s/regularization/result_%s_q%d-%d_lcurve_%s_reg_%.0e', p.optimization_output_path, p.sample_name, projection(1).par.r_sum{1}(1),  projection(1).par.r_sum{1}(end),...
        p.add_name, p.regularization_angle_coeff);
    print(fig_save, filename,'-dpng','-r300');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,~,proj_SH,Ereg] = optimization.SAXS_tomo_3D_err_metric([],p, s, projection(1:p.skip_projections:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save the calculated error
    E = 0;
    for ii = 1:numel(proj_SH)
        E = E + proj_SH(ii).error;
    end
    errormetric = [errormetric,E];
    penalty_error  =[penalty_error,Ereg];
    regularization_coeff = [regularization_coeff,p.regularization_angle_coeff];
end
toc(ticnow);

% plot L-curve
figure(40);
subplot(1, 2, 1)
hold on
plot(errormetric,(penalty_error./regularization_coeff),'ko--')
xlabel('error')
ylabel('regularization error')
ytickformat('%.2f')
xtickformat('%.2f')
% plot stabilization parameter (Glatter)
fig = subplot(1, 2, 2);
%hold on
left_color = [0 0 0];
right_color = [1 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
semilogx(regularization_coeff,penalty_error./regularization_coeff,'ko--')
grid on
xlabel('regularization parameter \mu')
ylabel('regularization error')
ytickformat('%.2f')
hold on
yyaxis right
semilogx(regularization_coeff,errormetric,'ro--')
grid on
ylabel('error')
ytickformat('%.2f')

% save
filename = sprintf('%s/regularization/find_mu_%s_q%d-%d_%s.mat', p.optimization_output_path, p.sample_name, projection(1).par.r_sum{1}(1),  projection(1).par.r_sum{1}(end),...
    p.add_name);
fprintf('Saving to: %s \n', filename)
save(filename,'errormetric', 'penalty_error', 'regularization_coeff', 'p', '-v6');

filename = sprintf('%s/regularization/find_mu_%s_q%d-%d_%s', p.optimization_output_path, p.sample_name, projection(1).par.r_sum{1}(1),  projection(1).par.r_sum{1}(end),...
    p.add_name);
fprintf('Saving figure to: %s \n', filename)
print(gcf, filename,'-dpng','-r300');