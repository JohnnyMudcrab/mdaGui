%% Function photons = computeNumOfPhotons2D (im, bgNoise, CCDSensitivity, EMGain)%, QE)
%
% computes the number of photons that contributed to the measurement of a gaussian peak depending on
% the camera settings
% NOTE: the computation doesn't consider the varius sources of noise
%
% parameters:
% -im				- the gaussian peak intensities
% -bgNoise			- background noise (measured in values of the image)
% -CCDSensitivity	- A/D converter discretization step size
% -EMGain			- EMGain used
%
% returns:
% -photons			- number of photons contributing to the peak
function photons = computeNumOfPhotons2D (im, bgNoise, CCDSensitivity, EMGain)%, QE)

% remove background noise from the intensity


% !!!!!!!!! THINK ABOUT THIS: i should not remove the background because this would set the lambda
% value of my poisson distributed intensities to a wrong value
% NOTE: if one does not want to remove the background simply provide a zero value for that parameter
im = im-bgNoise;
% !!!!!!!!!

% set all intensities to zero that are below zero now
im(im<0) = 0;

% compute number of photo electrons counted by the A/D converter
im = im.*CCDSensitivity;

% compute number of photo electrons before the EMGain thus those recognized by the ccd chip
im = im./EMGain;

% NOTE: we don't include the QE because those photons not creating a photo electron do not
% contribute to a measurement
% % compute number of photons according to the quantum efficiency of the camera
% im = im./QE;

% sum number of photons
photons = sum(sum(im));

end