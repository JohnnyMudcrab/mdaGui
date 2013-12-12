%% function fwhm1D = biasedFwhm1DDiscrete (pConfx, pConfy, x, dv)
%
% computes the difference function between two 1D fwhm dependences on z (parabola), which is required to compute
% the calibration curve from the different fwhm created by a cylindrical lens but can handle nonsteady grids which means that
% x,y pairs do not have to be provided linearly
% NOTE: after schuetz2000
%
% parameters:
% -pConfx		- vector of parameters defining the initial guess of the parabola in x direction:
%  		   -(1) - fwhm in focus (thus the resolution AND the angular point of the parabol)
%          -(2) - depth of focus
%		   -(3) - shifting of the parabol at the z-scale
%		   -(4) - gradient (m) of the straight line (mx+n) that biases the fwhm (no offset n)
%		   -..	- further parameters are up to you but are ignored here
% -pConfy		- vector of parameters defining the initial guess of the parabola in y direction
% -x      - x values (axial z)
% -dv		- real dimensions of one voxel in z direction
% -doCombinedCLensFit - if 1 the fwhm shift has to be considered else the curves are clearly defined and it is ignored
%
% returns:
% -fwhm1DDiffFcn		- y values of the difference function of two parabola (fwhm off focus)
%				NOTE: the fwhm will be returned in its real entity


% NOTE: that the difference function is only UNAMBIGUOUS if the two fwhm curves cross only ONCE and this only happens if they have the same configuration and
% there are ONLY shifts!!!!!!!!!!!!!!!!!!!!!!!!!!!
function fwhm1DDiffFcn = biasedFwhmDifferenceFcn1DDiscrete (pConfx, pConfy, x, dv, doCombinedCLensFit)

%%% compute x part

% create grids for the exponential calculation of every coordiante without using loops
% compute x vector
numeratorX = (x-pConfx(3)).*dv;

% compute the biased fwhm1D
fwhm1Dx = (pConfx(1) .* sqrt(1 + (numeratorX./pConfx(2)).^2 )) + pConfx(4).*numeratorX;

% in case this function is used with mda() and combined clens fit was used we add the fwhm offset
% if length(pConfx) > 5
if doCombinedCLensFit == 1
	fwhm1Dx = fwhm1Dx-pConfx(6);
end

%%% compute y part

% create grids for the exponential calculation of every coordiante without using loops
% compute y vector
numeratorY = (x-pConfy(3)).*dv;

% compute the biased fwhm1D
fwhm1Dy = (pConfy(1) .* sqrt(1 + (numeratorY./pConfy(2)).^2 )) + pConfy(4).*numeratorY;

% compute difference function
fwhm1DDiffFcn = fwhm1Dx-fwhm1Dy;

% figure,
% hold on
% plot(fwhm1Dx, 'b')
% plot(fwhm1Dy, 'r')
% plot(fwhm1Dx-fwhm1Dy, 'g')
% hold off
end