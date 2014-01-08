%% function [x1, x2] = inverseBiasedFwhm1D (pConf, y, dv)
%
% computes a 1D fwhm dependence on z (parabola) but can handle nonsteady grids which means that x,y pairs do not have to be
% provided linearly
% NOTE: after schuetz2000
%
% parameters:
% -pConf		- vector of parameters defining the initial guess of the parabola:
%  		   -(1) - fwhm in focus (thus the resolution AND the angular point of the parabol)
%          -(2) - depth of focus
%		   -(3) - shifting of the parabola at the z-scale
%		   -(4) - gradient (m) of the straight line (mx+n) that biases the fwhm (no offset n)
%		   -..	- further parameters are up to you but are ignored here
% -y		- the y vector (vector because than you can compute more results at once)
% -dv		- real axial dimension of one voxel
%
% returns:
% x		-	the x values corresponding to the y values where always two x for one y input is given
%			and the results are sorted in the second dimension for every y



function x = inverseBiasedFwhm1D (pConf, y, dv)

% check if there could be a division by zero but so far just warn
if isequal(pConf(4)^2 - (pConf(1)^2)/(pConf(2)^2), 0)
	fprintf('WARNING: division by zero in inverseBiasedFwhm1D\n');
end

% compute the according x values
% NOTE: because this parabol is unsymetric we realy have an x vector with 2 results
x(:,1) = ( (2.*y.*pConf(4) + sqrt( (2.*y.*pConf(4)).^2 - 4.*( pConf(4)^2 - (pConf(1)^2)/(pConf(2)^2) ) .* (y.^2-pConf(1)^2) ) )...
	./ ( 2* ( pConf(4)^2 - (pConf(1)^2)/(pConf(2)^2) ) ) ) + pConf(3)*dv;

x(:,2) = ( (2.*y.*pConf(4) - sqrt( (2.*y.*pConf(4)).^2 - 4.*( pConf(4)^2 - (pConf(1)^2)/(pConf(2)^2) ) .* (y.^2-pConf(1)^2) ) )...
	./ ( 2* ( pConf(4)^2 - (pConf(1)^2)/(pConf(2)^2) ) ) ) + pConf(3)*dv;

% check if there could be a negative number under the sqrt
if (sum(((2.*y.*pConf(4)).^2 - 4.*( pConf(4)^2 - (pConf(1)^2)/(pConf(2)^2) ) .* (y.^2-pConf(1)^2)<0))>0)

	% than set the according y values to the minima of the parabola which means that the value under
	% the root cannot become negative anymore
	
	% compute the x position of the minima of the parabola
	minX = sqrt( (pConf(4)^2*pConf(2)^4) / (pConf(1)^2- (pConf(4)^2*pConf(2)^2) ));
	% correct the direction
	if pConf(4)>0
		minX = -minX;
	end
	% add the axial shift
 	minX = minX+dv*pConf(3);
	% compute minima of the parabola
	minY = biasedFwhm1DDiscrete(pConf, minX./dv, dv);
	
	% recompute the according x values with the new y and set the value under root to zero
	x(((2.*y.*pConf(4)).^2 - 4.*( pConf(4)^2 - (pConf(1)^2)/(pConf(2)^2) ) .* (y.^2-pConf(1)^2)<0),1) ...
		= ( (2.*minY.*pConf(4) )...
		./ ( 2* ( pConf(4)^2 - (pConf(1)^2)/(pConf(2)^2) ) ) ) + pConf(3)*dv;
	x(((2.*y.*pConf(4)).^2 - 4.*( pConf(4)^2 - (pConf(1)^2)/(pConf(2)^2) ) .* (y.^2-pConf(1)^2)<0),2) ...
		= ( (2.*minY.*pConf(4) )...
		./ ( 2* ( pConf(4)^2 - (pConf(1)^2)/(pConf(2)^2) ) ) ) + pConf(3)*dv;
	
% 	if (2.*minY.*pConf(4)).^2 - 4.*( pConf(4)^2 - (pConf(1)^2)/(pConf(2)^2) ) .* (minY.^2-pConf(1)^2) <0
% 		fprintf('still negative :(\n');
% 		(2.*minY.*pConf(4)).^2 - 4.*( pConf(4)^2 - (pConf(1)^2)/(pConf(2)^2) ) .* (minY.^2-pConf(1)^2)
% 	end

	%%% NOTE: the following uncommented code proves that the value under the root can only become
	%%% negative if y is lower than the minima of the parabola and thus we can in this case just set
	%%% the root-value to zero
% 	fprintf('WARNING: negative sqrt in inverseBiasedFwhm1D\n');
%  	null = (2.*y.*pConf(4)).^2 - 4.*( pConf(4)^2 - (pConf(1)^2)/(pConf(2)^2) ) .* (y.^2-pConf(1)^2);
% 	null = null<0;
% 	null = y(null);
% 	
% 	nullX = sqrt( (pConf(4)^2*pConf(2)^4) / (pConf(1)^2- (pConf(4)^2*pConf(2)^2) ));
% 	if pConf(4)>0
% 		nullX = -nullX;
% 	end
%  	nullX = nullX+dv*pConf(3);
% 	nullY = biasedFwhm1DDiscrete(pConf, nullX./dv, dv);
% 	
% 	for i = 1:length(null)
% 		% test if the current y value causing the negative root is a value smaller than the smallest
% 		% possible value of the parabola
% 		if null(i)>nullY
% 			% if not
% 			fprintf('shit, the root can be negative even when the current y value is not smaller than the minima\n');
% 			null(i)
% 			nullY
% 		end
% 	end
end

% sort the x values
x = sort(x, 2);
% if x2<x1
% 	null = x2;
% 	x2 = x1;
% 	x1 = null;
% end

% %%% create grids for the exponential calculation of every coordiante without using loops
% % compute x and y vectors
% numeratorX = (x-pConf(3)).*dv;
% 
% % compute the biased fwhm1D
% fwhm1D = (pConf(1) .* sqrt(1 + (numeratorX./pConf(2)).^2 )) + pConf(4).*numeratorX;

end