function dens_reconstruct_all

    addpath('/mnt/das-gpfs/home/ext-georgiadis_m/Matias/scripts/utilities'); %path for utility functions

    path_data='/mnt/das-gpfs/home/ext-georgiadis_m/Matias/reconstruction'; % path to read and save datas and results for reconstruction
    load([path_data '/angles.mat']); % lookup table for rotation angles and scannumber, format is N*3 matrix of (rotx, roty, scannumber)
    load([path_data '/proj_all.mat']);
    load([path_data '/corr_mask.mat']);
    
    
    numiter=10000000;                                          % number of iterative corrections to run
    ratio=0.01;                                             % ratio of error correction, approximately 0.01, optimal value can be selected from testing with error list
    tomofilename=[path_data '/tomotrans_all.mat'];              % name of output model file
    modelsize=[427 427 495];
    subset=[300:598,600:880,882:1140];                                              % if to reconstruct from only a subset of all projections, leave empty to use all projections

    
    numprojs=size(angles,1);   
    if exist(tomofilename,'file')==0
        tomotrans=zeros(modelsize);
    else
        load(tomofilename);
    end
    if isempty(subset)
        subset=[1:numprojs];
    end
    
    for iter=1:numiter
        i=randi(size(subset,2));
        i=subset(i);
        rotx=angles(i,1);
        roty=angles(i,2);
        prawt=proj_all{i};
        psimt=proj_simul(tomotrans,rotx,roty,size(prawt));
        diff=(prawt-psimt).*corr_mask{i}.*ratio;
        tomotrans=proj_dens_correct(tomotrans,rotx,roty,diff,size(prawt));
        if mod(iter,100)==0||iter==1
            pshow=proj_simul(tomotrans,0,0,[495 427]);
            figure(1);
            imagesc(pshow);
            caxis([0 150]);
            pause(0.05);
            save(tomofilename,'tomotrans');
        end
        disp(iter);
    end
    save(tomofilename,'tomotrans');
end