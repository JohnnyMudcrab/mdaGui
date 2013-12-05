% function photons= computeNumOfParticlePhotons2D (im, dv, pos, fwhm, theta, bgNoiseMean, CCDSensitivity, EMGain)
%
% computes the number of photons that make up the peaks in the image
%
% parameters:
% -im			- the image containing the particles
% -dv			- real dimensions of one voxel
% -pos			- vector of xy positions of the center of the particles
% -fwhm			- vector of fwhms in xy of the PSFs fitted to the peaks
% -theta		- angle for the gaussian if an elliptical one is used (put [] for non-elliptical)
% -photonNumInterval - interval in number of sigmas for the gaussian that is used to collect the number of photons of the peak
% -CCDSensitivity	- A/D converter discretization step size
% -EMGain			- EMGain used
%
% returns (in nm!!!):
% -photons		- number of photons contributing to each peak
%
function [photons gIntensities] = computeNumOfParticlePhotons2D (im, dv, pos, fwhm, theta, photonNumInterval, bgNoiseMean, CCDSensitivity, EMGain)

%% parameters/variables

% % interval in number of sigmas for the gaussian that is used to collect the number of photons of the
% % peak
% photonNumInterval = 1.5;

%% preliminaries

%% compute number of photons for every peak
photons = [];

for i=1:size(pos, 1)
	%%% get number of photons for every peak
	
	% compute intensities overlapped by the peak
	if ~isempty(theta)
		gIntensities = getParticleIntensities2D(im, [1;1], pos(i, :), fwhm(i, :), theta(i), photonNumInterval);
	else
		gIntensities = getParticleIntensities2D(im, [1;1], pos(i, :), fwhm(i, :), theta, photonNumInterval);
	end
	
	% compute number of photons contributing to the peak
	tmpPhotons = computeNumOfPhotons2D(gIntensities, bgNoiseMean, CCDSensitivity, EMGain);

	% create list of number of photons
	photons = [photons; tmpPhotons];
	
end

end