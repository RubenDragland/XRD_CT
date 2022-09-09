function [cfg, vectors] = ...
                ASTRA_initialize(Npix, Nlayers,width_sinogram,angles,lamino_angle, tilt_angle, pixel_scale)
    %% make configs needed for astra, all angles are expected in degress 
            
    assert(math.isint(Npix), 'Npix is not integer');
    assert(math.isint(Npix), 'Nlayers is not integer');
    
    %% angles are assumed in degrees 
    if nargin < 5
        lamino_angle = 90; % tilt of the rotation axis in the direction of the beam (lamonigraphy angle)
    end
    if nargin < 6
        tilt_angle = 0; % rotation in the detector plane 
    end
    if nargin < 7
        pixel_scale = 1;  % px size on the detector [horizontal , vertical]
    end
    source_distance =1;
    show =  false; 
        
    Nangles = length(angles); 
    cfg.iVolX = Npix; 
    cfg.iVolY = Npix; 
    cfg.iVolZ = Nlayers(1); 
    cfg.iProjAngles = Nangles; 
    cfg.iProjU = width_sinogram; 
    cfg.iProjV = Nlayers(min(2,end)); % Nlayers can be a vector, useful for laminigraphy 
    cfg.iRaysPerDet = 1; 
    cfg.iRaysPerDetDim = 1; 
    cfg.iRaysPerVoxelDim = 1; 
    cfg.lamino_angle = lamino_angle;
    vectors = ASTRA_get_geometry(angles, lamino_angle, tilt_angle, source_distance,pixel_scale,show); 
    %%%% apply geometry correction to shift reconstruction into center %%%%%%%%%%%%%%%%%%%%%%%%%%
    vectors(:,4:6) = vectors(:,4:6) -(vectors(:,10:12)*(cfg.iProjV/2)+vectors(:,7:9)*(cfg.iProjU/2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     matlab2astra =  @(x)reshape(permute(x, [2,3,1]),[cfg.iProjU,cfg.iProjV,cfg.iProjAngles]); 
%     astra2matlab = @(x)permute(reshape(x, [cfg.iProjU,cfg.iProjAngles,cfg.iProjV]), [3,1,2]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function vectors = ASTRA_get_geometry(angles, lamino_angle, tilt_angle, source_distance,pixel_scale,show)

    Nangles = numel(angles);
    angles = deg2rad(angles);
    lamino_angle = pi/2 - deg2rad(lamino_angle); 
    tilt_angle  = deg2rad(tilt_angle); 
    if isscalar(lamino_angle)
        lamino_angle = lamino_angle .* ones(Nangles,1); 
    end
    if isscalar(tilt_angle)
        tilt_angle = tilt_angle .* ones(Nangles,2); 
    end
    if isscalar(pixel_scale) || numel(pixel_scale) == 2
        pixel_scale = bsxfun(@times, pixel_scale , ones(Nangles,2)); 
    end
        
    % We generate the same geometry as the circular one above. 
    vectors = zeros(Nangles, 12);
    for i = 1:Nangles
          % ray direction
          vectors(i,1) = sin(angles(i))*cos(lamino_angle(i));
          vectors(i,2) = -cos(angles(i))*cos(lamino_angle(i));
          vectors(i,3) = sin(lamino_angle(i));
                    
          vectors(i,1:3) =  vectors(i,1:3) *source_distance;
          % center of detector
          vectors(i,4:6) = 0;
          % vector from detector pixel (0,0) to (0,1)  
          vectors(i,7) = cos(angles(i))/pixel_scale(i,1);
          vectors(i,8) = sin(angles(i))/pixel_scale(i,1);
          vectors(i,9) = 0/pixel_scale(i,1);

          % vector from detector pixel (0,0) to (1,0)  
          % cross(vectors(i,1:3), vectors(i,7:9))
          % dot(vectors(i,1:3), vectors(i,7:9))

          vectors(i,10) = - sin(lamino_angle(i))*sin(angles(i))/pixel_scale(i,2);
          vectors(i,11) =   sin(lamino_angle(i))*cos(angles(i))/pixel_scale(i,2);
          vectors(i,12) =   cos(lamino_angle(i))/pixel_scale(i,2);     

        %  Rodrigues' rotation formula - rotate detector in plane
        %  perpendicular to the beam axis 
          vectors(i,7:9)=vectors(i,7:9)*cos(tilt_angle(i)) + ...
                         cross(vectors(i,1:3), vectors(i,7:9))*sin(tilt_angle(i)) + ...
                        (vectors(i,1:3)*dot(vectors(i,1:3),vectors(i,7:9)))*(1-cos(tilt_angle(i)));
          vectors(i,10:12)=vectors(i,10:12)*cos(tilt_angle(i)) + ...
                         cross(vectors(i,1:3), vectors(i,10:12))*sin(tilt_angle(i)) + ...
                        (vectors(i,1:3)*dot(vectors(i,1:3),vectors(i,10:12)))*(1-cos(tilt_angle(i)));
                   
                    
          
%       %% PLOT THE CURRENT SETUP 

      if show 
            clf
            draw_projection_geometry(vectors(i,:))
            pause(0.1);
      end



    end

end

function draw_projection_geometry(vectors)
    % show geometry saved in the "vectors" matrix 
    ray = vectors(1:3);
    c_center = vectors(4:6)-ray;
    c_origin = c_center - vectors( 7:9)/2-vectors( 10:12)/2;

    k = 6;
    n = 2^k-1;
    [x,y,z] = sphere(n);
    c = hadamard(2^k);
    s = 0.5;
    surf(vectors(4)+s*x,vectors(5)+s*y,vectors(6)+s*z,c);
    shading flat 
    colormap([1  1  0; 0  1  1])

    hold all
    plot3d_vec(c_origin, vectors( 7:9), 'r');
    plot3d_vec(c_origin, vectors( 10:12), 'r');
    plot3d_vec(c_origin+vectors( 7:9), vectors( 10:12), 'r');
    plot3d_vec(c_origin+vectors( 10:12), vectors( 7:9), 'r');
    % draw "pixels" on a 10x10 grid 
    for x = linspace(0,1,10)
        plot3d_vec(c_origin+x*vectors( 10:12), vectors( 7:9), 'r:');
        plot3d_vec(c_origin+x*vectors( 7:9),vectors( 10:12), 'r:');
    end
    mArrow3(c_center+ray*2,c_center,  'color', 'blue', 'stemWidth',0.02,'facealpha',0.5);

    hold off 
    axis([-1,1,-1,1,-1,1])
    
end
