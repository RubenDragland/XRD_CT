% saxs_tomo.m

addpath ~/Data10/cxs_software/base/

base_path='~/Data10/';
plot_sinograms=1;
savedata=1;

%for tomo
angular_intervals=360;
samx_intervals=400;
samy_intervals=3;
first_scan_no=4878;
last_scan=first_scan_no+(samy_intervals+1)*(angular_intervals+1)-1;

angle_step=180/angular_intervals;
angles=[0:angle_step:180];
angles_spec=zeros(angular_intervals+1,samy_intervals+1);
scans=[first_scan_no:last_scan];

first_scans=[first_scan_no:angular_intervals+1:last_scan];
last_scans=first_scans+angular_intervals;
no_of_line_intervals=zeros(size(first_scans))+samy_intervals;

%% Read integrated data and select a certain q region

stack=zeros(samy_intervals+1,samx_intervals,angular_intervals+1);
figure(1);
counter=1;
for jj=1:samy_intervals+1
    diode_mcs_filename=sprintf('%sanalysis/online/stxm/data/stxm_scans_%05d-%05d_mcs.mat',base_path,first_scans(jj),last_scans(jj));
    load(diode_mcs_filename);
    stack(jj,:,:)=mcs_data(2,1,:,:);
    if plot_sinograms
        imagesc(squeeze(stack(jj,:,:))); axis xy;
        title(sprintf('sinogram slice %03d',jj));
        drawnow;pause(0.2)
    end
end

% %% Plot movie with projections
% figure(1)
% for ii = 1:angular_intervals+1;
%     imagesc(stack(:,:,ii),[0,1e4]); axis xy equal tight;
%     title(sprintf('projection %03d',ii));
%     drawnow;pause(0.2);
% end

%% Plot sinogram and tomogram
% take offset from a corner of the sinogram:
length=10;
air_region=sum(sum(sum(stack(:,1:length,1:length))))/length/length/(samy_intervals+1);
tomo_filter='Hann';
freq_cutoff=1;
%beta=-log(stack-offset+1); % watch factor!
%beta=-log(stack);
beta=-log(stack/air_region);
figure2_position=[588   419   704   420];
figure(2);
set(gcf,'Position',figure2_position);
figure3_position=[20   269   560   570];
figure(3);
set(gcf,'Position',figure3_position);

for row=1:samy_intervals+1
    this_first_scan=first_scans(row);
    this_last_scan=last_scans(row)-1;
    title_str=sprintf('S%05d to S%05d',this_first_scan,this_last_scan);
    sino=squeeze(beta(row,:,1:end-1));
    sino_original=squeeze(stack(row,:,1:end-1));
    figure(2)
    imagesc(abs(sino)); %caxis([0 1e4])
    xlabel('angle'); ylabel('samx')
    title(sprintf('sinogram %s slice %02d',title_str,row));
    figure(3)
    slice = iradon(sino,angles(1:end-1), 'spline',tomo_filter,freq_cutoff,floor(size(sino,1)));
    imagesc(abs(slice)); axis xy equal tight
    title(sprintf('tomogram %s slice %02d',title_str,row));
    %caxis([0 300])
    pause(0.5)
    savefilename=sprintf('%sanalysis/saxs_tomo_diode/S%05d-%05d_slice_%03d',base_path,this_first_scan,this_last_scan,row);
    if savedata
        save([savefilename '.mat'],'sino','slice','sino_original','this_first_scan','this_last_scan','tomo_filter','freq_cutoff')
        print('-f3','-djpeg','-r300',[ savefilename '.jpg'] );
        print('-f3','-depsc','-r1200',[savefilename '.eps'] );
    end
end

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
%
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing la this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS scanning SAXS package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% Additionally, any publication using the package, or any translation of the 
%     code into another computing language should cite:
%    O. Bunk, M. Bech, T. H. Jensen, R. Feidenhans'l, T. Binderup, A. Menzel 
%    and F Pfeiffer, “Multimodal x-ray scatter imaging,” New J. Phys. 11,
%    123016 (2009). (doi: 10.1088/1367-2630/11/12/123016)
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
