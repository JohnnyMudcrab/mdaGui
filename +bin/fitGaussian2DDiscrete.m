%% Function gConf = fitGaussian2DDiscrete (im, dv, gConfInit, gConfLowerBound, gConfUpperBound, gType, weightedLsqPoisson)
%
% fits a gaussian function to a 2D image
% also allows different x and y width but can handle nonsteady grids which means
% that x,y,value triples do not have to be provided linearly
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
% -gConf		- vector representing the z values of the 2D gaussian at the positions given by the x
% and y vectors
function [gConf residual]= fitGaussian2DDiscrete (im, dv, gConfInit, gConfLowerBound, gConfUpperBound, gType, weightedLsqPoisson)

% just for communication with gaussian2DFrontend()
global voxelDim;
voxelDim = dv;

% just for communication with gaussian2DFrontend()
global exchangegType;
exchangegType = gType;

% just for communication with gaussian2DFrontend()
global exchangeWeightedLsqPoisson;
exchangeWeightedLsqPoisson = weightedLsqPoisson;

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

% sort out NaN values that are not used for the fitting
% create 2D mesh
numeratorX = 1:size(im, 2);
numeratorY = (1:size(im, 1))';
meshX=numeratorX(ones(length(numeratorY), 1) ,:);
meshY=numeratorY(:, ones(length(numeratorX), 1));

% remove nan values and create vectors from 2D meshes
x = (meshX(~isnan(im)))';
y = (meshY(~isnan(im)))';

% create z value vectors of gaussian from the image data
z = (im(~isnan(im)))';

% just for communication with gaussian2DFrontend()
global exchangeZ;
exchangeZ = z;

% start fitting routine

% if no weighting
if isempty(weightedLsqPoisson)
	[gConf null residual] = lsqcurvefit(@gaussian2DFrontend, gConfInit, [x;y], z...
		, gConfLowerBound, gConfUpperBound, optimset('Display', 'off'));
% weigth according to a poisson distribution
else
	[gConf null residual] = lsqcurvefit(@gaussian2DFrontend, gConfInit, [x;y], z.*sqrt(1./sqrt(z))...
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
% frontend to gaussian2DDiscrete() to be able use it with lsqcurvefit()
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
% -xy      - vector containing the xy positions at which the gaussian shall be computed
%
% returns:
% -gauss2D		- vector representing the z values of the 2D gaussian at the positions given by the x
% and y vectors
function gauss2D = gaussian2DFrontend (gConf, xy)

global voxelDim;
global exchangegType;
global exchangeWeightedLsqPoisson;
global exchangeZ;

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
		gauss2D = gaussian2DEllipseDiscrete(gConf, voxelDim, xy(1, :), xy(2, :));
	else
		gauss2D = gaussian2DDiscrete(gConf, voxelDim, xy(1, :), xy(2, :));
	end
% weigth according to a poisson distribution
else
	% if elliptical gaussian
	if exchangegType == 2
		gauss2D = gaussian2DEllipseDiscrete(gConf, voxelDim, xy(1, :), xy(2, :)).*sqrt(1./sqrt(exchangeZ));
	else
		gauss2D = gaussian2DDiscrete(gConf, voxelDim, xy(1, :), xy(2, :)).*sqrt(1./sqrt(exchangeZ));

% 		% update the intensity value of gConf because for holtzer07 this is required as the number of photons and would be too low otherwise
% 		gConf(5) = (gConf(5)*pi*sqrt((gConf(3)*voxelDim(1))^2*(gConf(4)*voxelDim(2))^2)) / (4*log(2));
% 		gauss2D = aberratedGaussian2DDiscrete(gConf, voxelDim, xy(1, :), xy(2, :)).*sqrt(1./sqrt(exchangeZ));

	end
end

end