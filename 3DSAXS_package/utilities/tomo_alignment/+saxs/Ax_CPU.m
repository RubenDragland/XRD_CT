function projData = Ax_CPU(volData, cfg, vectors)


    %% default parameters, FIXME 
    p.volume_upsampling = 1;
    p.filter_2D = 3;           % Lower values give you compromise between resolution and artifacts
    p.method = 'bilinear';     % 'nearest' or 'bilinear'

    Npix = cfg.iVolX; 
    grid = -Npix/2+1 :Npix/2;
    zgrid =  -cfg.iVolZ/2+1 :cfg.iVolZ/2;
    [X,Y,Z] = meshgrid(grid, grid, zgrid);


    volData = double(volData);

    tilt_angle = 3/2*pi-atan2(vectors(:,3), sqrt(sum(vectors(:,1:2).^2,2))); 
    theta = 3*pi/2 -atan2(vectors(:,1)./cos(tilt_angle), vectors(:,2) ./ cos(tilt_angle)); 
    
    ugrid = -cfg.iProjU/2+1:cfg.iProjU/2;
    vgrid = -cfg.iProjV/2+1:cfg.iProjV/2;

    
    projData = zeros(cfg.iProjV, cfg.iProjU, cfg.iProjAngles);
    parfor ii = 1:cfg.iProjAngles
        R = saxs.rotation_matrix_3D(  tilt_angle(ii), 0,   theta(ii));
        projData(:,:,ii) = arb_projection( volData, X, Y, Z, R, p, ugrid, vgrid );
    end
    


end