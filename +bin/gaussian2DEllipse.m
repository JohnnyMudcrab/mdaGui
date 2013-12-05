%% Function gauss2D = gaussian2DEllipse (gConf, dv, sizex, sizey)
%
% computes a 2D elliptical gaussian that can be turned but need a full 2D grid as input
% NOTE: this is the most GENERAL form of a gaussian function!!!!!!!!!!!!!!!
%
% parameters:
% -gConf	- vector of parameters defining the gauss:
%      -(1) - center x position
%      -(2) - center y position
%      -(3) - fwhm x
%      -(4) - fwhm y
%      -(5) - height of the peak
%      -(6) - offset
%      .. reserved
%	   -(9) - turning angle theta if an elliptical gaussian is used
%	   -..	- further parameters are up to you but are ignored here
% -dv		- xyz dimension of one voxel
% -sizex	- size of the underlying grid in x direction
% -sizey	- size of the underlying grid in y direction
%
% returns:
% gauss2D	- matrix representing the 2D gaussian

function gauss2D = gaussian2DEllipse (gConf, dv, sizex, sizey)

% convert fwhm into a parameter scaling the width of the gaussian function using:
% sigma = fwhm/ 2*sqrt(2*log(2))
sigmaX = gConf(3)*dv(1) / 2.3548200450309;
sigmaY = gConf(4)*dv(2) / 2.3548200450309;

% convert angle theta to radian
theta = gConf(9)*pi/180;

% compute matrix parameters
a = ((cos(theta)^2)/(2*sigmaX^2)) + ((sin(theta)^2)/(2*sigmaY^2));

b = -(sin(2*theta)/(4*sigmaX^2)) + (sin(2*theta)/(4*sigmaY^2));

c = ((sin(theta)^2)/(2*sigmaX^2)) + ((cos(theta)^2)/(2*sigmaY^2));

%% vectorized solution

%%% create grids for the exponential calculation of every coordiante without using loops
% compute x and y vectors
numeratorX = (dv(1):dv(1):(sizex*dv(1)))-(gConf(1)*dv(1));
numeratorY = (dv(2):dv(2):(sizey*dv(2)))'-(gConf(2)*dv(2));
% create matrixes to be able to vectorize the problem (shift the values into the other rows/columns)
meshX=numeratorX(ones(length(numeratorY), 1) ,:);
meshY=numeratorY(:, ones(length(numeratorX), 1));

% compute the gaussian 2D
gauss2D = gConf(6) + gConf(5) .* exp( -( (a.*meshX.^2) + (2.*b.*meshX.*meshY) + (c.*meshY.^2) ) );

end