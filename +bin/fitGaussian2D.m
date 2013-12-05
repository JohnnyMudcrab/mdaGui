%% Function gConf = fitGaussian2D (im, dv, gConfInit, gConfLowerBound, gConfUpperBound, gType, weightedLsqPoisson)
%
% fits a gaussian function to a 2D image
%
% parameters:
% -im			- the image on which the gauss shall be fitted to
% -dv			- real dimensions of one voxel
% -gConfInit	- vector of parameters defining the initial guess of the gauss:
%        -(1)	- center x position
%        -(2)	- center y position
%        -(3)	- fwhm x
%        -(4)	- fwhm y
%        -(5)	- height of the peak
%        -(6)	- offset
%        .. reserved
%	     -(9) - turning angle theta if an elliptical gaussian is used
%	     -..	- further parameters are up to you but are ignored here
% -gConfLowerBound	- vector of parameters defining the lower bound of the configuration of the	gauss
% -gConfUpperBound	- vector of parameters defining the upper bound of the configuration of the gauss
% -gType		- (0) -> gaussian with possibly varying width in x and y direction
%				- (1) -> gaussian with the same width in x and y direction
%				- (2) -> elliptical gaussian having varying x and y widths but can also be turned by
%				an angle theta
% -weightedLsqPoisson	- if defined than a weighting parameter adjusted for poisson distributions
% is used for the lsq agorithm
%
% returns:
% -gConf		- vector of parameters of the fitted gauss

function [gConf residual] = fitGaussian2D (im, dv, gConfInit, gConfLowerBound, gConfUpperBound, gType, weightedLsqPoisson)

% just for communication with gaussian2DFrontend()
global voxelDim;
voxelDim = dv;

% just for communication with gaussian2DFrontend()
global exchangegType;
exchangegType = gType;

% just for communication with gaussian2DFrontend()
global exchangeWeightedLsqPoisson;
exchangeWeightedLsqPoisson = weightedLsqPoisson;

% just for communication with gaussian2DFrontend()
global exchangeIm;
exchangeIm = im;

% if the width in xy shall remain the same
if gType == 1
	% delete the y fwhm parameter from the configurations
	gConfInit(4) = [];
	gConfLowerBound(4) = [];
	gConfUpperBound(4) = [];
end

% if elliptical gaussian
if gType == 2
	% save and delete the absolute xy locations (7,8) from gConfs
	savedgConfInit = gConfInit(7:8);
	gConfInit(7:8) = [];
	gConfLowerBound(7:8) = [];
	gConfUpperBound(7:8) = [];
end

% if no weighting
if isempty(weightedLsqPoisson)
	[gConf null residual] = lsqcurvefit(@gaussian2DFrontend, gConfInit, size(im), im...
		, gConfLowerBound, gConfUpperBound, optimset('Display', 'off'));
% weigth according to a poisson distribution
else
	[gConf null residual] = lsqcurvefit(@gaussian2DFrontend, gConfInit, size(im), im.*sqrt(1./sqrt(im))...
		, gConfLowerBound, gConfUpperBound, optimset('Display', 'off'));
end

% if the width in xy shall remain the same
if gType == 1
	% finally add the same fwhm for y as we have fitted for x
	gConf = horzcat(gConf(1:3), gConf(3), gConf(4:end));
end

% if elliptical gaussian
if gType == 2
	% insert the old values for the absolute xy locations (7,8) again
	gConf = horzcat(gConf(1:6), savedgConfInit, gConf(7:end));
end


end

%% Function gaussian2DFrontend (gConf, sizexy)
%
% frontend to gaussian2D() to be able use it with lsqcurvefit()
%
% parameters:
% -gConf		- vector of parameters defining the initial guess of the gauss:
%  		   -(1) - center x position
%		   -(2) - center y position
%          -(3) - fwhm x
%          -(4) - fwhm y NOTE: depending on exchangegType
%          -(5) - height of the peak
%          -(6) - offset
%		   -..	- further parameters are up to you but are ignored here
% -sizeyx      - vector containing the xy size of the underlying image
%
% returns:
% -gauss2D		- vector of parameters of the fitted gauss

function gauss2D = gaussian2DFrontend (gConf, sizeyx)

global voxelDim;
global exchangegType;
global exchangeWeightedLsqPoisson;
global exchangeIm;

% if the width in xy shall remain the same
if exchangegType == 1
	% just add the x fwhm for the y fwhm in gConf as well such that we don't have to change the
	% gaussian function
	gConf = horzcat(gConf(1:3), gConf(3), gConf(4:end));
end

% if elliptical gaussian
if exchangegType == 2
	% just add the absolute xy positions again such that we don't have to change the
	% gaussian function
	gConf = horzcat(gConf(1:6), [0 0], gConf(7:end));
end

% if no weighting
if isempty(exchangeWeightedLsqPoisson)
	% if elliptical gaussian
	if exchangegType == 2
		gauss2D = gaussian2DEllipse(gConf, voxelDim, sizeyx(2), sizeyx(1));
	else
		gauss2D = gaussian2D(gConf, voxelDim, sizeyx(2), sizeyx(1));
	end
% weigth according to a poisson distribution
else
	% if elliptical gaussian
	if exchangegType == 2
		gauss2D = gaussian2DEllipse(gConf, voxelDim, sizeyx(2), sizeyx(1)).*sqrt(1./sqrt(exchangeIm));
	else
		gauss2D = gaussian2D(gConf, voxelDim, sizeyx(2), sizeyx(1)).*sqrt(1./sqrt(exchangeIm));
	end
end

end