function [gConfs] = locateInvoker(img, threshold, pixelSizeX, pixelSizeY, FWHMmin, FWHMmax, mask)

  trackingParams{1}  = 10;
  trackingParams{2}  = 0;% 0.5; % 4; %synapseDistInterval
  trackingParams{3}  = [FWHMmin,...
                        FWHMmin;
                        FWHMmax,...
                        FWHMmax];
  trackingParams{4}  = 0;% gType
  trackingParams{5}  = 4;
  trackingParams{6}  = 25;% minNumPhotonsPerPeak
  trackingParams{7}  = 4;% fitQualityThresholdFactor
  trackingParams{8}  = 2/3;% gaussianPart2Fit
  trackingParams{9}  = 0.1;% intensityOutliers
  trackingParams{10} = [0 1 0 0 1];% reasonabilityChecks
  trackingParams{11} = 40;% 3; %1;% subImSizefactor
  trackingParams{12} = false;
  trackingParams{13} = 0;% use weighted LSQ or not
  trackingParams{14} = [300 300];% [400 400] initFwhm
  trackingParams{15} = 1.5;% watershedRegionControlAxisSizeFactor
  trackingParams{16} = 1;% gaussianOffsetConstraints, 1: std-bounded, 2: fixed or 3: unbounded gaussian offset
  trackingParams{17} = 0;%0; %1;% removeAllParticlesTooCloseToOthers 0: no 1: yes
  % setting: useImTopHatFilter = [iDiffRel, dTop, dBrim] == [relative intensity difference (as a factor for the std), diameter of top and brim]
  trackingParams{18} = [4 3 7];%[4 3 7];%[0.5 9 13];%[];% useImTopHatFilter
  trackingParams{19} = 0; %- doEMFitting -> 0 NOTE: if set to >0 (100) than the useWatershedTransformRegions flag is ignored and treated as 0 thus a certain rectangular region is
  %		used for fitting
  trackingParams{20} = 1; % useConfForSPT -> if set to 0 then all is as usual and the definition of parameters makes sense for detecting synapses with varying sizes
  %						if set to 1 then the code is optimized for SPT where all particles are basically similar and where the definition of distances
  %						and the size of the fitting regions are set constant and do not depend on the size of a particle anymore.
  %						thus the parameters: 2 (synapseDistInterval) is now interpreted as the absolute distance in pixels that two particles must have
  trackingParams{21} = 1; % doSimpleButRobustThresholding -> if set to 1 than the new thresholding using the median is used and also the median is used for the imTopHatFilter
  trackingParams{22} = 1; % useHighResolutionTopHatFilter -> if set to 1 than the new topHatFiler is used that is more robust (by using the mean instead of max) and allows to distinguish nearby synapses (NOTE: its slower)					and 11 (subImSizeFactor) is an absolute value in pixels for the size of the fitting region


  [gConfs, ~, ~] = locateParticles2DBetha (img, [pixelSizeX pixelSizeY], trackingParams, ...
                                           0, 0, [], [], mask, []);


end

