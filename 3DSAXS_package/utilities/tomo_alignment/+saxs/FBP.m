function [rec,sinogram] = FBP(sinogram, cfg, vectors, split, valid_angles)

    if nargin < 4
        split = 1;
    end

    if nargin > 4 && ~isempty(valid_angles)
        sinogram = sinogram(:,:,valid_angles);
        vectors = vectors(valid_angles,:);
        try cfg.lamino_angle = cfg.lamino_angle(valid_angles); end 
    end
        
    
    [Nlayers,Nw,Nproj] = size(sinogram);
%     assert(cfg.iProjAngles == Nproj, 'Wrong number of angles')
    cfg.iProjAngles = Nproj;
    assert(cfg.iProjU == Nw, 'Wrong sinogram width')
    assert(cfg.iProjV == Nlayers, 'Wrong sinogram height')
    assert(mod(Nw,1)==0, 'Only even width of sinogram is supported')


    split_fft = prod(split);
    split_fft= max(split_fft, ceil(numel(sinogram)*4/1e9)); % split to 1GB blocks 
    if isa(sinogram, 'gpuArray') && split_fft>1
        gpu = gpuDevice; 
        fprintf('Free GPU memory: %3.2g%%\n',  gpu.AvailableMemory/gpu.TotalMemory*100)       
        sinogram = fft_partial(sinogram,2,1,split_fft);
    else
        sinogram = fft(sinogram,[],2);
    end
        
    % !!! call the variable by the same name to avoid memory copy 
    %% MATLAB FILTER 
    win = 1;  % hann(Nw)'
    d = 1;
    filter = 'ram-lak';
    filter = designFilter(filter, Nw, d);
    filter = bsxfun(@times, filter', reshape(sind(cfg.lamino_angle),1,1,[]));
    sinogram = bsxfun(@times, sinogram, filter);
     if isa(sinogram, 'gpuArray')  && split_fft>1
        sinogram = real(ifft_partial(sinogram,2,1, split_fft));
    else
        sinogram = real(ifft(sinogram,[],2));
    end
    sinogram = sinogram * (pi/2/Nproj);  % constant to make the scale identical with the input phantom !!! 
    
    
    %% unfiltered backprojection 
%     try
%         parallel.gpu.GPUDevice.isAvailable
%         use_gpu = true; 
%     catch
%         use_gpu = false; 
%     end
        
    use_gpu = false; 

    if use_gpu
        rec = tomo.Atx_partial(sinogram, cfg, vectors,split);
    else       
        rec = saxs.Atx_CPU(sinogram, cfg, vectors);
    end

    
    
    
end



%======================================================================
function filt = designFilter(filter, len, d)
% Returns the Fourier Transform of the filter which will be
% used to filter the projections
%
% INPUT ARGS:   filter - either the string specifying the filter
%               len    - the length of the projections
%               d      - the fraction of frequencies below the nyquist
%                        which we want to pass
%
% OUTPUT ARGS:  filt   - the filter to use on the projections


% order = max(64,2^nextpow2(2*len));
order = len;

if strcmpi(filter, 'none')
    filt = ones(1, order);
    return;
end

% First create a bandlimited ramp filter (Eqn. 61 Chapter 3, Kak and
% Slaney) - go up to the next highest power of 2.

n = 0:(order/2); % 'order' is always even. 
filtImpResp = zeros(1,(order/2)+1); % 'filtImpResp' is the bandlimited ramp's impulse response (values for even n are 0)
filtImpResp(1) = 1/4; % Set the DC term 
filtImpResp(2:2:end) = -1./((pi*n(2:2:end)).^2); % Set the values for odd n
filtImpResp = [filtImpResp filtImpResp(end-1:-1:2)]; 
filt = 2*real(fft(filtImpResp)); 
filt = filt(1:(order/2)+1);

% keyboard

w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist

switch filter
    case 'ram-lak'
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
    case 'cosine'
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
    case 'hamming'
        filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
    case 'hann'
        filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
    otherwise
        error(message('images:iradon:invalidFilter'))
end

filt(w>pi*d) = 0;                      % Crop the frequency response
filt = [filt' ; filt(end-1:-1:2)'];    % Symmetry of the filter
%----------------------------------------------------------------------
end
