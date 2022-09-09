function volData = Atx_CPU(projData, cfg, vectors)


    %% default parameters, FIXME 
    p.volume_upsampling = 0;
    p.filter_2D = 0;           % Lower values give you compromise between resolution and artifacts
    p.method = 'bilinear';     % 'nearest' or 'bilinear'

    Npix = cfg.iVolX; 
    grid = -Npix/2+1 :Npix/2;
    zgrid =  -cfg.iVolZ/2+1 :cfg.iVolZ/2;
    [X,Y,Z] = meshgrid(grid, grid, zgrid);


    projData = double(projData);

%     projData = padarray(projData,[0,(cfg.iVolX-cfg.iProjU)/2,0],0,'both'); 
    
    tilt_angle = 3/2*pi-atan2(vectors(:,3), sqrt(sum(vectors(:,1:2).^2,2))); 
    theta = 3*pi/2 -atan2(vectors(:,1)./cos(tilt_angle), vectors(:,2) ./ cos(tilt_angle)); 
    
    ugrid = -cfg.iProjU/2+1:cfg.iProjU/2;
    vgrid = -cfg.iProjV/2+1:cfg.iProjV/2;
        
    volData = zeros(cfg.iVolX,cfg.iVolY,cfg.iVolZ);
    parfor ii = 1:cfg.iProjAngles
        R = (saxs.rotation_matrix_3D(  tilt_angle(ii), 0,   theta(ii)));
        volData = volData + arb_back_projection( projData(:,:,ii), ugrid,vgrid, X, Y, Z, R, p );
    end



end