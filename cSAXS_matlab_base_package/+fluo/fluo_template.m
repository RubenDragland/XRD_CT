%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example template for fluorescence analysis functionalities %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot spectra of scan without calibration
 
scans = 735;
ii = 1;
falcon = io.nexus_read('~/Data10','scanNr',scans(ii),'filter','falcon');

falcon_spectra = squeeze(fluo.fluo_raw_to_spectra(falcon));
aux = load('+fluo/detector_calibration_200915_163627.mat');
aux.calibration

figure(3)
plot(log10(mean(falcon_spectra,2)));
% plot(log10(falcon_spectra));
title(sprintf('Scan %s',scans))

