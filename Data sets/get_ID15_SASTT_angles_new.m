function [rot_x, rot_y] = get_ID15_SASTT_angles_new()


fname = 'SASTT_phantom4_19Feb2020.mat';
f = load(fname);

for ii = 1:length(f.projection)
    
    rot_x(ii) = f.projection(ii).rot_x;
    rot_y(ii) = f.projection(ii).rot_y;
end
    

end