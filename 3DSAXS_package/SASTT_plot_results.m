%% load aligned projection file and the result from optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_path = '/das/work/p17/p17633/adrian/'; % path for offline analysis
sample_name = 'tomo_1';                     %'sample_name';% name given in the saxs_caller_template
add_name = 'std';%'ID';       % additional name the optimizations: = [ ] if not needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath utilities/
filename = fullfile(base_path,sprintf('analysis/SASTT/%s/projection_data/SASTT_%s_aligned_FBP.mat', ...
    sample_name, sample_name));
orig_data = load(filename);
% load the optimization results from SH
filename = fullfile(base_path, sprintf('analysis/SASTT/%s/SH/%s/results/', sample_name, add_name));
        filename = fullfile(filename, sprintf('result_%s_q%d-%d_all_again_%s.mat', sample_name, ...
            orig_data.projection(1).par.r_sum{1}(1), ...
            orig_data.projection(1).par.r_sum{1}(end), add_name));
load(filename);
%% plot view3
tomo = s.a(1).data;             %3D array
% viewpar                       Structure with other parameters for visualization
viewpar.coloraxis = 'auto'      %Axis for the colorbar, e.g. = [-1 1] or 'auto';default 'auto'
viewpar.which = [];             %Vector with position where cuts are displayed, default [round(sizey/2) round(sx/2) round(sz/2)]
viewpar.interactive = false;     %If (= true, default) you can click the slices and
                                %see the corresponding perpendicular slices corresponding to the point you
                                %clicked. If (= false) then it just plots once.
%%%
view3(tomo,viewpar)
%% plot view3 with spherical harmonics in each point

% viewpar                       Structure with other parameters for visualization
viewpar.coloraxis = 'auto'      %Axis for the colorbar, e.g. = [-1 1] or 'auto';default 'auto'
viewpar.which = [];             %Vector with position where cuts are displayed, default [round(sizey/2) round(sx/2) round(sz/2)]
viewpar.interactive = true;     %If (= true, default) you can click the slices and
                                %see the corresponding perpendicular slices corresponding to the point you
                                %clicked. If (= false) then it just plots once.
viewpar.which_coeff =2;         %which coeffiecient to plot, default = 1
viewpar.plotradius = 0;         %If false the function is shown on a sphere, if true the
                                %radius is proportional to the function
%%%%%%
view3_SASTT(filename,viewpar)
%% plot spherical harmonics squared from specific point in tomogram (as it is in the error metric of SASTT)
%select y,x,z from tomogram to plot as spherical harmonics (eg using view3
x = 23;
y = 11;
z = 24;
plotradius = 0;                 %for 1 scale the radius of spherical harmonics with its intensity
figure_nummer = 50;
%%%%%%%%%
figure(figure_nummer)
plot_spherical_harmonic_squared([s.theta.data(y,x,z),s.phi.data(y,x,z)],...
    horzcat(s.a.l),horzcat(s.a.m),...
    [s.a(1).data(y,x,z),s.a(2).data(y,x,z),s.a(3).data(y,x,z),s.a(4).data(y,x,z)],...
    plotradius)
colormap jet
%% plot spherical harmonics with own values
orientation =  [pi/2,pi/2];     %[theta phi], Rotation of axes is performed first by theta
                                %along y and then by phi along z. in radians
l = [0 2 4];                    %Vector of polar orders
m = [0 0 0];                    %Vector of azimuthal orders
a = [1 0.5 0.2];               %Vector of coefficients, same length as l and m
plotradius = 0;                 %If false the function is shown on a sphere, if true the
                                %radius is proportional to the function
squared = true;                 %If true plot square of spherical harmonics,
                                %as it is in error metric
figure_nummer = 51;
%%%%%%%%%%
figure(figure_nummer)
if squared
    plot_spherical_harmonic_squared(orientation,l,m,a,plotradius)
    colormap jet
else
    plot_spherical_harmonic(orientation,l,m,a,plotradius)
    colormap jet
end
%% plot 3D angle
slicenumber = [];               %(y,x,z) Number of slice to display if empty plot middle
%%%
plot_3D_color(s,slicenumber)

%% look at transmission and scattering from all projections
figure_no = 70;
for ii = 1:length(orig_data.projection)
    figure(figure_no)
    
    imagesc(orig_data.projection(ii).diode)
    axis xy
    axis equal tight
    colormap bone
    colorbar
    title(sprintf('transmission %d',ii))
   
    figure(figure_no+1)
    
    imagesc(sum(orig_data.projection(ii).data,3))
    axis xy
    axis equal tight
    colormap bone
    colorbar
    title(sprintf('scattering %d',ii))
    drawnow
    pause(0.1)
    end
%% look at azimuthal integration from one point
x = 23;
y = 11;
projection_nr = 1;
figure_no = 80;
%%%%%%%%%%%%%
figure(figure_no)
plot(orig_data.projection(projection_nr).integ.phi_det(1:8),squeeze(orig_data.projection(projection_nr).data(y,x,:)))
xlabel('azimuthal angle [deg]')
ylabel('Intensity')