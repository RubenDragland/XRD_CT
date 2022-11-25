function R = calcRotExp(rot_x,rot_y,fliprot_x,fliprot_y)

    if fliprot_x
        rot_x = -rot_x;
    end
    if fliprot_y
        rot_y = -rot_y;
    end

%     %%%%% ORIGINAL VERSION
%     %calculate rotation matrix
%     %Rotation of object around x: tilt axis A
    rotangle = deg2rad(rot_x); % transform to radians
    R =     [1      0               0;
        0  cos(rotangle)   -sin(rotangle);
        0  sin(rotangle)   cos(rotangle)];
    %Rotation of object around y: tomo axis B
    rotangle = deg2rad(rot_y); % transform to radians
    R =R * [cos(rotangle) 0  sin(rotangle);
        0      1        0      ;
        -sin(rotangle) 0  cos(rotangle)];
%     
    

% %%%% BASED ON ACTA CRYST AND WOLFRAM DEFINITIONS
% 
%     rotangle = deg2rad(rot_x); % transform to radians
%     R =     [1      0               0;
%         0  cos(rotangle)   sin(rotangle);
%         0  -sin(rotangle)   cos(rotangle)];
% 
%     rotangle = deg2rad(rot_y); % transform to radians
%     R =R * [cos(rotangle) 0  -sin(rotangle);
%         0      1        0      ;
%         sin(rotangle) 0  cos(rotangle)];
%     
    
    
    
end

