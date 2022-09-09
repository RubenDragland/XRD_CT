%% Loading data
addpath ..
addpath ../ptycho/utils/
addpath ../ptycho/

clear

% close all;

scans =  [91:146]; % scan numbers
fign = 100;
filename_list={};
scans_found=[];

for ii=1:numel(scans)
    try
        ptycho_filename = find_ptycho_filename('~/Data10/analysis/',scans(ii),[],'recons');
        if ~iscell(ptycho_filename)
            filename_list=[filename_list,utils.abspath(ptycho_filename)];
        else
            filename_list=[filename_list,utils.abspath(ptycho_filename{end})];
        end
        scans_found = [scans_found,scans(ii)];
    catch err
        disp(err)
    end

end
    scans = scans_found;
    p0=io.HDF.hdf5_load(filename_list{1},'/reconstruction/p');
    
    Nthreads = 6;
    
    dims_ob=p0.object_size;
    pixsize=p0.dx_spec(1);
    size_crop=[round((p0.object_size(1)-p0.asize(1))/2)*2,round((p0.object_size(2)-round(p0.asize(2)))/2)*2];

    
    proj_file_names  = reshape(filename_list, 1,[]); 
    
    object_block = io.ptycho_read(Nthreads, 'single', dims_ob, '/reconstruction/object', proj_file_names); 
    object_block=permute(object_block,[2 1 3]);
    size_now=size(object_block);
    pixsize=zeros(numel(scans),1);
    
    for ii=1:numel(scans)
        energy=io.HDF.hdf5_load(filename_list{ii},'/reconstruction/p/energy');
        energies(ii)=energy;
        dx_spec=io.HDF.hdf5_load(filename_list{ii},'/reconstruction/p/dx_spec');
        pixsize(ii)=dx_spec(1);
    end
%         object_resize=utils.imrescale_fft(object_block(:,:,ii),dx_spec(1)/pixsize);
    for ii=1:numel(scans)
        object_resize=utils.imrescale_fft(object_block(round(size_now(1)/2)-size_crop(1)/2:round(size_now(1)/2)+size_crop(1)/2,round(size_now(2)/2)-size_crop(2)/2:round(size_now(2)/2)+size_crop(2)/2,ii),pixsize(ii)/max(pixsize));
        object_crop(:,:,ii)=object_resize;
    end
    
amp=abs(object_crop);
phase=math.unwrap2D_fft2(object_crop,10,0,[],-1);


%% Selecting the nonair part
roi_chosen = {1:323,15:431};   % (y,x)
% If the sample is in the middle of FOV, include the entire object so only
% air is outside this ROI.
% If air region is in the middle of FOV, include only air.

obj_or_air=1; % =1 if object is in middle, =0 if air is in middle.

figure(fign+5); clf;
plotting.imagesc3D(phase); 
colormap bone; axis xy equal tight;
caxis([-2*pi 0.2]);
hold on;
plotting.vline(roi_chosen{2}(1),'-b');
plotting.vline(roi_chosen{2}(end),'-b');
plotting.hline(roi_chosen{1}(1),'-b');
plotting.hline(roi_chosen{1}(end),'-b');

drawnow;

figure(fign+6); clf;
plotting.imagesc3D(amp); 
colormap bone; axis xy equal tight;
caxis(math.sp_quantile(amp,[0.01 0.99],10));
hold on;

plotting.vline(roi_chosen{2}(1),'-b');
plotting.vline(roi_chosen{2}(end),'-b');
plotting.hline(roi_chosen{1}(1),'-b');
plotting.hline(roi_chosen{1}(end),'-b');

size_t=size(amp);

if obj_or_air==1
    mask=ones(size_t(1),size_t(2));
    mask(roi_chosen{:})=0;
else
    mask=zeros(size_t(1),size_t(2));
    mask(roi_chosen{:})=1;
end

for ii=1:numel(scans)
    trans = abs(object_crop(:,:,ii));
    reg(ii)=mean(mean(trans(mask==1)));
    trans = trans/reg(ii);
    transmission (:,:,ii) = trans;
end

%% Alignment of projections
object_crop = utils.stabilize_phase(object_crop); 

phase_unwrap=math.unwrap2D_fft2(object_crop,50,0,[],-1); %[roi_chosen{2}(1) size(phase,2)-roi_chosen{2}(2)] );
% ref=mean(phase_unwrap,3);
ref=phase_unwrap(:,:,1);

shift = utils.find_shift_fast_2D(phase_unwrap, ref, 0.01); 
% shift = utils.dftregistration(fft2(ref),fft2(phase_unwrap(:,:,2)),100); 
phase_unwrap = utils.imshift_fft(phase_unwrap, -shift); 
transmission_aligned = utils.imshift_fft(transmission, -shift);



figure(fign+5); clf; 
plotting.imagesc3D(phase_unwrap); 
colormap bone; axis xy equal tight;
caxis([-2*pi 0.2]);
hold on;
plotting.vline(roi_chosen{2}(1),'-b');
plotting.vline(roi_chosen{2}(end),'-b');
plotting.hline(roi_chosen{1}(1),'-b');
plotting.hline(roi_chosen{1}(end),'-b');
drawnow;

figure(fign+6); clf; 

plotting.imagesc3D(transmission_aligned); 
colormap bone; axis xy equal tight;
caxis(math.sp_quantile(transmission_aligned,[0.01 0.99],10));
hold on;
plotting.vline(roi_chosen{2}(1),'-b');
plotting.vline(roi_chosen{2}(end),'-b');
plotting.hline(roi_chosen{1}(1),'-b');
plotting.hline(roi_chosen{1}(end),'-b');

%% XANES analysis

roi = {90:190;100:290};  %(y,x)

figure(fign+5); clf; 
plotting.imagesc3D(phase_unwrap,'fps',2); 
colormap bone; axis xy equal tight;
caxis([-2*pi 0.2]);
hold on;
plotting.vline(roi{2}(1),'-r');
plotting.vline(roi{2}(end),'-r');
plotting.hline(roi{1}(1),'-r');
plotting.hline(roi{1}(end),'-r');

plotting.vline(roi_chosen{2}(1),'-b');
plotting.vline(roi_chosen{2}(end),'-b');
plotting.hline(roi_chosen{1}(1),'-b');
plotting.hline(roi_chosen{1}(end),'-b');

figure(fign+fign+6); clf; 
plotting.imagesc3D(transmission_aligned,'fps',2); 
colormap bone; axis xy equal tight;
caxis(math.sp_quantile(transmission_aligned,[0.01 0.99],10));
hold on;
plotting.vline(roi{2}(1),'-r');
plotting.vline(roi{2}(end),'-r');
plotting.hline(roi{1}(1),'-r');
plotting.hline(roi{1}(end),'-r');

plotting.vline(roi_chosen{2}(1),'-b');
plotting.vline(roi_chosen{2}(end),'-b');
plotting.hline(roi_chosen{1}(1),'-b');
plotting.hline(roi_chosen{1}(end),'-b');


if ~isempty(roi)
    xanes_phase = squeeze(mean(mean(phase_unwrap(roi{:},:))));
    xanes_amp   = squeeze(mean(mean(transmission_aligned(roi{:},:))));
else
    xanes_phase = squeeze(mean(mean(phase_unwrap)));
    xanes_amp   = squeeze(mean(mean(transmission_aligned)));
end
figure(fign+7); plot(energies,-xanes_phase(:),'ro-'); %.*energies(:).^2);
title('phase')
figure(fign+8); plot(energies,xanes_amp(:),'ro-');
title('amp')
figure(fign+9); plot(energies(1:end-1),diff(xanes_amp(:)));
title('amp\_diff')

figure(fign+10); 
plot(energies,shift*max(pixsize)*1e6);
legend('\Deltax','\Deltay')
xlabel('E [keV]')
ylabel('Displacement [\mum]')

%% EXAFS

[r,m,b]=regression(energies(40:end),xanes_amp(40:end)');
exafs=1./(xanes_amp(30:end))-1./((m*energies(30:end)'+b));
ex_k = sqrt((energies(30:end)-6.557)*1000*0.2625);
k=[0:0.1:max(ex_k)];

exafs_k=interp1(ex_k,(exafs-0.2*exp(-ex_k')).*(ex_k').^2,k);
%%

dr=pi/numel(exafs_k)/0.1; 
exafs_k(isnan(exafs_k))=0;
figure; plot([1:numel(exafs_k)]*dr,abs(fft(exafs_k)'));