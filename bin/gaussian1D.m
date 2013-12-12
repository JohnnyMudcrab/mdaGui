%% Function gauss1D = gaussian1D (gConf, dv, sizex)
%
% computes a 1D gaussian
%
% parameters:
% -gConf	- vector of parameters defining the gauss:
%      -(1) - center x position
%      -(2) - fwhm
%      -(3) - height of the peak
%      -(4) - offset
%	   -..	- further parameters are up to you but are ignored here
% -dv		- real dimensions of one voxel
% -sizex	- size of the underlying grid in x direction
%
% returns:
% gauss1D	- matrix representing the 1D gaussian

function gauss1D = gaussian1D (gConf, dv, sizex)

% convert fwhm into a parameter scaling the width of the gaussian function using:
% sigma = fwhm/ 2*sqrt(2*log(2))
sigma = gConf(2)*dv / 2.3548200450309;

%%% create grids for the exponential calculation of every coordiante without using loops
% compute x and y vectors
numeratorX = ((dv:dv:(sizex*dv))-(gConf(1)*dv)).^2;
% % create matrixes to be able to vectorize the problem (shift the values into the other rows/columns)
% meshX=numeratorX(ones(length(numeratorX), 1) ,:);

% compute the gaussian 1D
gauss1D = gConf(4) + gConf(3) .* exp( -numeratorX ./ (2*sigma^2));

end