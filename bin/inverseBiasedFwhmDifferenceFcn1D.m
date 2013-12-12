%% x = inverseBiasedFwhmDifferenceFcn1D(y, lookUpTable, dv)
%
% computes the inverse difference function between two 1D fwhm dependences on z (parabola), which is required to compute
% the calibration curve from the different fwhm created by a cylindrical lens but can handle nonsteady grids which means that
% x,y pairs do not have to be provided linearly
% NOTE: after schuetz2000
%
% parameters:
% -y		- y value you wnat the x position for
% -xGrid	- x grid for the lookup table
% -lookUpTable	-all y values for the x values used to get the belonging x value
%
% returns:
% -fwhm1DDiffFcn		- y values of the difference function of two parabola (fwhm off focus)
%				NOTE: the fwhm will be returned in its real entity

% NOTE: that the difference function is only UNAMBIGUOUS if the two fwhm curves cross only ONCE and this only happens if they have the same configuration and
% there are ONLY shifts!!!!!!!!!!!!!!!!!!!!!!!!!!!
function x = inverseBiasedFwhmDifferenceFcn1D(y, xGrid, lookUpTable)
	% get absolute x index of inverse function
	[null idx] = min(abs(lookUpTable-y));
	
	% compute index according to the underlying matrix grid
	x = xGrid(idx);
end