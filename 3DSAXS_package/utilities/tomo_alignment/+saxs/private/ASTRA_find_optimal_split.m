function split = ASTRA_find_optimal_split(cfg)
%     reset(gpuDevice)
    gpu = gpuDevice;
    
    Nangles = min(1024, cfg.iProjAngles);
    if isfield(cfg, 'Grouping')
        Nangles = min(Nangles, cfg.Grouping);
    end
    split = 1; 
    split = max(split, ceil([cfg.iVolX, cfg.iVolY, cfg.iVolZ] /  2000)); % texture memory limit
    split = max(split, [1,1,ceil( (cfg.iVolX*cfg.iVolY*cfg.iVolZ*4) /  1.024e9 / prod(split))]); % texture memory limit
    split = max(split, [1,1,ceil( (cfg.iVolX*cfg.iVolY*cfg.iVolZ*4*2) / gpu.AvailableMemory / prod(split))]); % gpu memory limit 
    split(4)  = ceil( (cfg.iProjU*cfg.iProjV*Nangles*8) / gpu.AvailableMemory); % gpu memory limit 
    
    split = 2.^nextpow2(split);
    if any(split ~= 1); fprintf('Automatically splitting to %ix%ix%ix(%i) cubes \n', split); end

end