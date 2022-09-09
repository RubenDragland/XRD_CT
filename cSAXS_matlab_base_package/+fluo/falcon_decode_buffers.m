%FALCON_DECODE_BUFFERS decode raw file produced by falcon
%
% ** buffer_data            raw data
%
% returns:
% ++ falcon_data            decoded dataset
%
% EXAMPLES:
%       % load raw h5 file and parse the buffer data into the decoder
%       fluo_data = io.HDF.hdf5_load('~/Data10/data/S00000-00999/S00250/e16812_1_00250.h5')
%       buffer_data = fluo_data.entry.data.data;
%       decoded_buffer = falcon_decode_buffers(buffer_data);
%
% IMPORTANT NOTE:
%       unfortunately, there is no information on the number of spectra
%       acquired in the map and the decoded buffer will have additional
%       spectra that are results of the way how falcon electronics
%       transfers the data (for example in chunks of 16 spectra). So it is
%       crucial to load positions data and extract number of points
%       expected to reduce the decoded dataset. This is automatically done
%       in the function fluo_read
%
% see also: fluo.fluo_read()

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a
%   publication or if it is fully or partially rewritten for another
%   computing language the authors and institution should be acknowledged
%   in written form in the publication: “Data processing was carried out
%   using the “cSAXS software package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.”
%   Variations on the latter text can be incorporated upon discussion with
%   the CXS group if needed to more specifically reflect the use of the package
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%
% This code and subroutines are part of a continuous development, they
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its
%    proper use and the correctness of the results.

function [falcon_data] = falcon_decode_buffers(buffer_data)
buffer_data = int16(buffer_data);  % the falcon data is saved in 16-bit words
                                   % and needs to be parsed in decoder for 
                                   % typecast function to work properly


%% default parameters
maxPixels = 0;
nPixels = 0;
clockPeriod = 320e-9;
firstTime = 1;
maxPixelsPerBuffer = 0;

%% prepare empty header structures
falcon_buffer_header.tag0 = [];  % Tag word 0
falcon_buffer_header.tag1 = [];  % Tag word 1
falcon_buffer_header.headerSize = []; % Buffer header size
falcon_buffer_header.mappingMode = []; % Mapping mode (1=Full spectrum, 2=Multiple ROI, 3=List mode)
falcon_buffer_header.runNumber = []; % Run number
falcon_buffer_header.bufferNumber = []; % Sequential buffer number, low word first
falcon_buffer_header.bufferID = []; % 0=A, 1=B
falcon_buffer_header.numPixels = []; % Number of pixels in buffer
falcon_buffer_header.startingPixel = []; % Starting pixel number, low word first
falcon_buffer_header.moduleNumber = [];
falcon_buffer_header.channelID = [];
falcon_buffer_header.channelElement = [];
falcon_buffer_header.reserved1 = zeros(6, 1);
falcon_buffer_header.channelSize = [];
falcon_buffer_header.reserved2 = zeros(3, 1);
falcon_buffer_header.bufferErrors = [];


falcon_mca_pixel_header.tag0 = []; % Tag word 0
falcon_mca_pixel_header.tag1 = []; % Tag word 1
falcon_mca_pixel_header.headerSize = []; % Buffer header size
falcon_mca_pixel_header.mappingMode = []; % Mapping mode (1=Full spectrum, 2=Multiple ROI, 3=List mode)
falcon_mca_pixel_header.pixelNumber = []; % Pixel number
falcon_mca_pixel_header.blockSize = []; % Total pixel block size, low word first
falcon_mca_pixel_header.spectrumSize = [];


falcon_data.firstPixel = [];
falcon_data.numPixels = [];
falcon_data.Data = [];
falcon_data.RealTime = [];
falcon_data.LiveTime = [];
falcon_data.InputCounts = [];
falcon_data.OutputCounts = [];


%% extract information on the buffer data
s = size(buffer_data);
nDimensions = numel(s);
bufferSize = s(1);

nDetectors = 1;
nArrays = 1;
if (nDimensions > 1), nDetectors = s(2); end
if (nDimensions > 2), nArrays = s(3); end


%% decoding of the falcon buffers
for array=1:nArrays
    for detector=1:nDetectors
        % Copy the buffer header information
        temp_data = buffer_data(:, detector, array);
        falcon_buffer_header.tag0 = temp_data(1);
        falcon_buffer_header.tag1 = temp_data(2);
        falcon_buffer_header.headerSize = temp_data(3);
        falcon_buffer_header.mappingMode = temp_data(4);
        falcon_buffer_header.runNumber = temp_data(5);
        falcon_buffer_header.bufferNumber = typecast(temp_data(6:7), 'int32') + 1;
        falcon_buffer_header.bufferID = temp_data(8);
        falcon_buffer_header.numPixels = temp_data(9);
        falcon_buffer_header.startingPixel = typecast(temp_data(10:11), 'int32') + 1;
        if array == 1, firstPixel = falcon_buffer_header.startingPixel; end
        falcon_buffer_header.moduleNumber   = temp_data(12);
        falcon_buffer_header.channelID      = temp_data(13);
        falcon_buffer_header.channelElement = temp_data(14);
        falcon_buffer_header.channelSize    = temp_data(21);
        falcon_buffer_header.bufferErrors   = temp_data(25);
        
        if falcon_buffer_header.mappingMode == 1 % MCA mapping mode
            offset = 256;
            % It is possible that each buffer has a different number
            % of pixels.  Use the value in the first buffer as the
            % maximum for all buffers
            if maxPixelsPerBuffer == 0, maxPixelsPerBuffer = falcon_buffer_header.numPixels; end
            if falcon_buffer_header.numPixels <= maxPixelsPerBuffer, nPix = falcon_buffer_header.numPixels; end
            if falcon_buffer_header.numPixels > maxPixelsPerBuffer, nPix = maxPixelsPerBuffer; end
            if detector == 1, nPixels = nPixels + nPix; end
            
            for pixel = 1:nPix
                pn = (array-1) * maxPixelsPerBuffer + pixel;
                falcon_mca_pixel_header.tag0        = temp_data(offset + 1);
                falcon_mca_pixel_header.tag1        = temp_data(offset + 2);
                falcon_mca_pixel_header.headerSize  = temp_data(offset + 3);
                falcon_mca_pixel_header.mappingMode = temp_data(offset + 4);
                falcon_mca_pixel_header.pixelNumber = typecast(temp_data(offset + 5:offset + 6), 'int32') + 1;
                % The following assumes that the netCDF file is one single
                % xMAP run. But that may not be true, so ignore the
                % pixelNumber value pn = mph.pixelNumber - firstPixel
                % Instead set it as above, which also handles case of
                % different number of pixels in each buffer
                
                falcon_mca_pixel_header.blockSize   = typecast(temp_data(offset + 7:offset + 8), 'int32');
                % The spectrum size appears to be in units of 16-bit words
                falcon_mca_pixel_header.spectrumSize = temp_data(offset + 9) / 2;
                if firstTime
                    firstTime = 0;
                    nChannels = falcon_mca_pixel_header.spectrumSize;
                    maxPixels = maxPixelsPerBuffer * nArrays;
                    data      = zeros(nChannels, nDetectors, maxPixels);
                    realTime     = zeros(nDetectors, maxPixels);
                    liveTime     = zeros(nDetectors, maxPixels);
                    inputCounts  = zeros(nDetectors, maxPixels);
                    outputCounts = zeros(nDetectors, maxPixels);
                end
                p = offset + 32;
                realTime(detector, pn)     = double(typecast(temp_data(p + 1:p + 2), 'int32')) * clockPeriod;
                liveTime(detector, pn)     = double(typecast(temp_data(p + 3:p + 4), 'int32')) * clockPeriod;
                inputCounts(detector, pn)  = typecast(temp_data(p + 5:p + 6), 'int32');
                outputCounts(detector, pn) = typecast(temp_data(p + 7:p + 8), 'int32');                
                p = offset + 256;
                channel_iterator_double = 1;
                for channel_iterator = 1 : nChannels
                    data(channel_iterator, detector, pn) = typecast(temp_data(p + ...
                         channel_iterator_double:p + 1 + channel_iterator_double), 'int32');
                    channel_iterator_double = channel_iterator_double + 2;
                end
                offset = offset + falcon_mca_pixel_header.blockSize;
            end % for pixel = 1:nPix
        end % if falcon_buffer_header.mappingMode == 1
    end  % for detector=1:nDetectors
end  % for array=1:nArrays

%% export data and information
if falcon_buffer_header.mappingMode ~= 3
    falcon_data.firstPixel   = firstPixel;
    falcon_data.numPixels    = nPixels;
    falcon_data.Data         = data;
    falcon_data.RealTime     = realTime;
    falcon_data.LiveTime     = liveTime;
    falcon_data.InputCounts  = inputCounts;
    falcon_data.OutputCounts = outputCounts;
end
end