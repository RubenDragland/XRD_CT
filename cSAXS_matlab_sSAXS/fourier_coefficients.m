function [f1_amp, f2_amp, f2_phase, I_cos_dev] = fourier_coefficients(data)

%if numel(size(data)) == 2
    num_segments = size(data, 1)*2; % multiplied by 2 because it is the number of segments on the detector
    % Fourier transform the intensity
    if size(num_segments > 1)
        I_fft = fft(data);
    else
        I_fft(1) = data;
        I_fft(2) = 0;
    end
%elseif numel(size(data)) == 3%
    %num_segments = size(data, 3)*2; % multiplied by 2 because it is the number of segments on the detector
    %I_fft = fft(data);
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Real instead of abs returns negative base lines with the
% correct sign. Negative base lines may occur in case of
% background subtraction.
% SYMMETRIC INTENSITY: f1
f1_amp =    real(I_fft(1)) ./(num_segments/2); % because we average based on symmetry: half detector
%ASYMMETRIC INTENSITY: f2
% The factor of two takes into account that the first and
% last Fourier component contribute equally.
f2_amp = 2 * abs(I_fft(2)) ./(num_segments/2); % because we average based on symmetry: half detector

% phase of f2
% The phase offset, i.e., orientation.% offset the zero-point of the rotation from the start of the first
% segment to its center
% From Oliver this correction was done in plot_fourier_result four_dat.f2_phase = four_dat.f2_phase - pi/8;
f2_phase = angle(I_fft(2));

%%%%%%%%%%%%%%%%%%%%
% Computing the RMS deviation using the residuum of the FFT
I_fft2 = I_fft;
I_fft2([1 2 end]) = 0;  % Setting to zero, zero and first coefficients
I_cos_dev = sqrt(mean( abs(  ifft(I_fft2)  ).^2 ));
