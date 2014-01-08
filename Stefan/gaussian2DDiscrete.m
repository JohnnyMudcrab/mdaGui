%% Function gauss2D = gaussian2DDiscrete (gConf, dv, x, y)
%
% computes a 2D gaussian that allows different x and y width but can handle nonsteady grids which means
% that x,y,value triples do not have to be provided linearly
%
% parameters:
% -gConf	- vector of parameters defining the gauss:
%      -(1) - center x position
%      -(2) - center y position
%      -(3) - fwhm x
%      -(4) - fwhm y
%      -(5) - height of the peak
%      -(6) - offset
%	   -..	- further parameters are up to you but are ignored here
% -dv		- xyz dimension of one voxel
% -x	- the x vector
% -y	- the y vector
%
% returns:
% gauss2D	- vector representing the z values of the 2D gaussian at the positions given by the x
% and y vectors

function gauss2D = gaussian2DDiscrete (gConf, dv, x, y)

% convert fwhm into a parameter scaling the width of the gaussian function using:
% sigma = fwhm/ 2*sqrt(2*log(2))
sigmaX = gConf(3)*dv(1) / 2.3548200450309;
sigmaY = gConf(4)*dv(2) / 2.3548200450309;

%% vectorized solution

%%% create grids for the exponential calculation of every coordiante without using loops
% compute x and y vectors
numeratorX = (x.*dv(1)-(gConf(1)*dv(1))).^2;
numeratorY = (y.*dv(2)-(gConf(2)*dv(2))).^2;
% % create matrixes to be able to vectorize the problem (shift the values into the other rows/columns)
% meshX=numeratorX(ones(length(numeratorX), 1) ,:);
% meshY=numeratorY(:, ones(length(numeratorY), 1));

% compute the gaussian 2D
gauss2D = gConf(6) + gConf(5) .* exp( (-numeratorX./(2*sigmaX^2)) - (numeratorY./(2*sigmaY^2)) );

% %% solution using for-loops
% 
% % allocate final matrix
% gauss2D = zeros(sizey, sizex);
% 
% for x=1:sizex
% 	for y=1:sizey
% 		gauss2D(y, x) = exp( -( ((x*dv(1))-gConf(1)*dv(1)).^2 )./ (2*sigmaX^2)...
% 			-( ((y*dv(2))-gConf(2)*dv(2)).^2 )./ (2*sigmaY^2) );
% 	end
% end
% 
% gauss2D = gConf(6) + gConf(5) .* gauss2D;

end