function dens_reconstruct(projection)

    
    numiter=1000;                                          % number of iterative corrections to run
    ratio=0.03;                                             % ratio of error correction, approximately 0.01, optimal value can be selected from testing with error list
    tomofilename=['~/Work/matlab/BT_tomogram.mat'];              % name of output model file
    projsize=size(projection(1).data);
    modelsize=[projsize(1),projsize(2),projsize(2)];
    subset=[];
    
    numprojs=length(projection);   

    if isempty(subset)
        subset=[1:numprojs];
    end
    
    proj_shifted=cell(numprojs);
    
    for ii = 1:length(projection)
        proj_shifted{ii} = double(utils.shiftwrapbilinear(mean(projection(ii).data,3).*projection(ii).window_mask, ...
            -projection(ii).dy, -projection(ii).dx));
    end
    
    tomotrans=zeros(modelsize);
    
    for iter=1:numiter
        i=randi(size(subset,2));
        i=subset(i);
        prawt=proj_shifted{i};
        psimt=proj_simul(tomotrans,projection(i).Rot_exp,size(prawt));
        diff=(prawt-psimt).*ratio;
        tomotrans=proj_dens_correct(tomotrans,projection(i).Rot_exp,diff,size(prawt));
        if 1 %mod(iter,100)==0||iter==1
            pshow=proj_simul(tomotrans,eye(3),size(prawt));
            figure(1);
            imagesc(pshow);
            drawnow;
%             save(tomofilename,'tomotrans');
        end
        disp(iter);
    end
    save(tomofilename,'tomotrans');
end