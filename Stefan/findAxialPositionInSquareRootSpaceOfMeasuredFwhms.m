% finds z (axial position) that has the minimal distance between the measured fwhms and the calibration curves in the square root space
% axialRange - holds the lower and upper bound of the axial range such that the zAxis can be computed (must be provided in nm)
% zAxisResolution - used to define the precision of the zAxis for finding closest point on curve (in nm)
function [z fwhmXCalib fwhmYCalib] = findAxialPositionInSquareRootSpaceOfMeasuredFwhms(zInit, fwhmX, fwhmY, fwhmConfX, fwhmConfY, dv, axialRange, zAxisResolution, doCombinedCLensFit)

global fwhmXGlobal;
global fwhmYGlobal;
global fwhmConfXGlobal;
global fwhmConfYGlobal;
global dvGlobal;
global doCombinedCLensFitGlobal;
global zAxisGlobal;
global fwhmXCurveGlobal;
global fwhmYCurveGlobal;
global zAxisResolutionGlobal;

fwhmXGlobal = fwhmX;
fwhmYGlobal = fwhmY;
fwhmConfXGlobal = fwhmConfX;
fwhmConfYGlobal = fwhmConfY;
dvGlobal = dv;
doCombinedCLensFitGlobal = doCombinedCLensFit;
zAxisResolutionGlobal = zAxisResolution;

% compute the calibration curves with high resolution (in nm)
% zAxisResolution = 0.1;
% zAxis = -527:zAxisResolution:300;
zAxisGlobal = axialRange(1):zAxisResolution:axialRange(2);
if doCombinedCLensFit
	fwhmXCurveGlobal = biasedFwhm1DDiscrete(fwhmConfX, zAxisGlobal./dv(3), dv(3))-fwhmConfX(6);
else
	fwhmXCurveGlobal = biasedFwhm1DDiscrete(fwhmConfX, zAxisGlobal./dv(3), dv(3));
end
fwhmYCurveGlobal = biasedFwhm1DDiscrete(fwhmConfY, zAxisGlobal./dv(3), dv(3));

% search for z that has the minimal distance between the measured fwhms and the calibration curves in the square root space
z = fminsearch(@squareRootMinDist, zInit);
if doCombinedCLensFit
	fwhmXCalib = biasedFwhm1DDiscrete(fwhmConfX, z/dv(3), dv(3))-fwhmConfX(6);
else
	fwhmXCalib = biasedFwhm1DDiscrete(fwhmConfX, z/dv(3), dv(3));
end
fwhmYCalib = biasedFwhm1DDiscrete(fwhmConfY, z/dv(3), dv(3));

end

%% subfunction that computes the distance of the fwhm values and the calibration curves in the square root space of the fwhms

function dist = squareRootMinDist(z)

global fwhmXGlobal;
global fwhmYGlobal;
global fwhmConfXGlobal;
global fwhmConfYGlobal;
global dvGlobal;
global doCombinedCLensFitGlobal;
global zAxisGlobal;
global fwhmXCurveGlobal;
global fwhmYCurveGlobal;
global zAxisResolutionGlobal;

%% similar to suppl. huang 2008 but computes the distance to the calibration curves as the euclidean distance

% discretize z according to the resolution of our z axis
z = round(z/zAxisResolutionGlobal)*zAxisResolutionGlobal;

% find the shortest euclidean distance between the current measurement and the calibration curves individually
distx = sqrt( (fwhmXGlobal-fwhmXCurveGlobal).^2 + (z-zAxisGlobal).^2 );
[null minxIdx] = min(distx);
disty = sqrt( (fwhmYGlobal-fwhmYCurveGlobal).^2 + (z-zAxisGlobal).^2 );
[null minyIdx] = min(disty);

dist = sqrt( ( distx(minxIdx) )^2 ...
	+		 ( disty(minyIdx) )^2 );
% z
% %% suppl. huang 2008 (computes the distance to the calibration curves only along the y axis)
% if doCombinedCLensFitGlobal
% 	dist = sqrt( ( sqrt(fwhmXGlobal) - sqrt(biasedFwhm1DDiscrete(fwhmConfXGlobal, z/dvGlobal(3), dvGlobal(3))-fwhmConfXGlobal(6)) )^2 ...
% 		+		 ( sqrt(fwhmYGlobal) - sqrt(biasedFwhm1DDiscrete(fwhmConfYGlobal, z/dvGlobal(3), dvGlobal(3))) )^2 );
% else
% 	dist = sqrt( ( sqrt(fwhmXGlobal) - sqrt(biasedFwhm1DDiscrete(fwhmConfXGlobal, z/dvGlobal(3), dvGlobal(3))) )^2 ...
% 		+		 ( sqrt(fwhmYGlobal) - sqrt(biasedFwhm1DDiscrete(fwhmConfYGlobal, z/dvGlobal(3), dvGlobal(3))) )^2 );
% end
% z
end