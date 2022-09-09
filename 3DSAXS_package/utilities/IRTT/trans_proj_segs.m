function trans_proj_segs(util,path_recon)
    addpath(util);
     load([path_recon , '/proj_segs_raw.mat']);
%     load([path_recon , '/delta_proj.mat']);
     load([path_recon , '/corr_mask.mat']);
    n=size(proj_segs_raw,2);
    proj_segs=cell(1,n);
    for i=1:n
        data=proj_segs_raw{i};
        out=zeros(size(data));
        for j=1:8
            out(j,:,:)=proj_translation(squeeze(data(j,:,:)).*corr_mask{i},delta_proj(i,:));
        end
        proj_segs{i}=out;
    end
    save([path_recon , '/proj_segs.mat'],'proj_segs');
end
