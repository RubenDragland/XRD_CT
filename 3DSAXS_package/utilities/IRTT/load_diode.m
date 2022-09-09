function load_diode
% Load diode signals that can be used to correct signal absorption in
% reconstruction.

    fileformat_diode='/mnt/das-gpfs/work/p16266/analysis/online/stxm/data/stxm_scans_%.5d-%.5d_mcs.mat';
    % Format of diode data output from SAXS scanning script.
    
    path_data='/mnt/das-gpfs/work/p16266/TOMO-SAXS/2_Reconstruction'; % path to read and save datas and results for reconstruction
    load([path_data '/angles.mat']); % lookup table for rotation angles and scannumber, format is N*3 matrix of (rotx, roty, scannumber)

    numprojs=size(angles,1);    % number of projections, should be the same with rotation angles
    numlines=70;                % number of lines (scans) in one projection
    numpoints=106;              % number of points in one line
    
    Norm_intensity=10000;       % a global scale to normalize the diode readout to values around 1
    
    
    % modify until here
    diodeabsorption=zeros(numprojs,numpoints,numlines);
    for i=1:numprojs
        firstscan=angles(i,3);
        load(sprintf(fileformat_diode,firstscan,firstscan+numlines-1));
        d=flip(squeeze(mcs_data(3,1,:,:)),1);
        if size(d,1)<numpoints
            for j=(size(d,1)+1):numpoints
                d(j,:)=d(j-1,:);
            end
        end
        for j=1:numlines
            if mod(j,2)==0
                d(:,j)=flip(d(:,j),1);
            end
        end
        d=d/Norm_intensity;
        diodeabsorption(i,:,:)=d(:,:);
    end
    save([path_data '/diode2.mat'],'diodeabsorption');
end
