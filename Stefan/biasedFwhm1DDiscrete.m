%% function fwhm1D = biasedFwhm1DDiscrete (pConf, x, dv)
%
% computes a 1D fwhm dependence on z (parabola) but can handle nonsteady grids which means that x,y pairs do not have to be
% provided linearly
% NOTE: after schuetz2000
%
% parameters:
% -pConf		- vector of parameters defining the initial guess of the parabola:
%  		   -(1) - fwhm in focus (thus the resolution AND the angular point of the parabol)
%          -(2) - depth of focus
%		   -(3) - shifting of the parabol at the z-scale
%		   -(4) - gradient (m) of the straight line (mx+n) that biases the fwhm (no offset n)
%		   -..	- further parameters are up to you but are ignored here
% -x      - x values (axial z)
% -dv		- real dimensions of one voxel
%
% returns:
% -fwhm1D		- y values of the parabola (fwhm off focus)
%				NOTE: the fwhm will be returned in its real entity


function fwhm1D = biasedFwhm1DDiscrete (pConf, x, dv)

%%% create grids for the exponential calculation of every coordiante without using loops
% compute x and y vectors
numeratorX = (x-pConf(3)).*dv;

% compute the biased fwhm1D
fwhm1D = (pConf(1) .* sqrt(1 + (numeratorX./pConf(2)).^2 )) + pConf(4).*numeratorX;

end