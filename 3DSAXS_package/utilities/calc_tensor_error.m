function err=calc_tensor_error(projection,tomotensor,B_segs)
    E=0;
    numOfpixels=length(projection)*size(projection(1).data,1)*size(projection(1).data,2);
    for ii=1:length(projection)
        for kk = 1:size(projection(ii).data, 3)
            proj_shifted(:,:,kk) = double(utils.shiftwrapbilinear(squeeze(projection(ii).data(:,:,kk)), ...
                -projection(ii).dy, -projection(ii).dx));
        end
        B=squeeze(B_segs(ii,:,:));
        projsize=size(proj_shifted);
        projsize=projsize([3,1,2]);
        proj_sim = proj_simul_tensor(tomotensor, projection(ii).Rot_exp, B, projsize);
        proj_sim(proj_sim<0)=0;
%         proj_sim = abs(proj_sim);
        proj_sim = permute(proj_sim,[2 3 1]);
        aux_diff_poisson = (sqrt(proj_sim) - sqrt(proj_shifted)).* utils.shiftwrapbilinear(projection(ii).window_mask, ...
                -projection(ii).dy, -projection(ii).dx);
        E = E + 2*sum(sum(sum(aux_diff_poisson.^2))) / numOfpixels;
    end
    err=E;
end