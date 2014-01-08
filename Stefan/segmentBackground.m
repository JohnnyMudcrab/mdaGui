% function FGMask = segmentBackground(im, bgWidth, debug)
%
% returns a binary mask with the background set to zero and the foreground
% set to one
% NOTE: the returned background standard deviation does not include the morphological changes made
% to close holes in the background. thus the returned value is smaller as if you would compute the
% standard deviation from all intensities that this function markes as background.

% -this function can handle 2D and 3D images
%
% parameters:
% -im			- the n-dimensional image
% -bgWidth		- factor, in terms of the standard deviation of the background distribution, that is used to define
%               which values still belong o the background
% -debug		- 1 for debugging
%
% returns:
% -FGMask		- binary mask with the background set to zero and the foreground set to one
% -bgMean		- mean value of the background
% -bgRMS		- standard deviation of the background (NOTE: this one might be too small)
function [FGMask bgMean bgRMS] = segmentBackground(im, bgWidth, debug)

%% parameters

% signal to background ratio, which is used to take only those range of the histogram into account
% for fitting the gaussian that contains the background
SBR = 2;

% offset for the range of the histogram to ignore the possibly huge amount of zero values that screw
% up the histogram fitting
histogramOffset = 2;

% factor for the peak index in the histogram, which is used to as a cutter for the histogram length
% to save computational costs as well as to be able to ignore signal outliers
histCutFactor = 3;

%% main

% compute the histgram
imVec = im(:);
% NOTE: +1 because the image can have zero values as well!
histogram = zeros(max(imVec)+1, 1);
for i = 1:length(imVec)
	histogram(imVec(i)+1) = histogram(imVec(i)+1) + 1;
end

% figure,plot(histogram);

% take only the lower half of the histogram range to ignore outlayers as
% well as to save computational costs
% NOTE: in case there is nearly no signal than either the histogram is not cutted or a factor from
% the maximal peak index is taken
% histogram = histogram(histogramOffset:round(length(histogram)/SBR));
[null maxIdx] = max(histogram);
histogram = histogram(histogramOffset:max(round(length(histogram)/SBR), min(length(histogram), maxIdx*histCutFactor)));

%% do a gaussfit on the histogram. this should find the distribution of the
% background pixels
% the gaussFit is computed by optimizing the least-squares difference
% between a gaussian curve and the histogram

% create construct for initial gaussian parametrization
% the mean is taken from the index in the histogram with the highest value
[null, meanInit] = max(histogram);
sigmaInit = 2.3548200450309*meanInit;
bgConfInit = [];
bgConfInit(1) = meanInit;
bgConfInit(2) = sigmaInit;
bgConfInit(3) = max(histogram);
bgConfInit(4) = 0;

bgConfInitLB = [];
bgConfInitLB(1) = 0;
bgConfInitLB(2) = 0;
bgConfInitLB(3) = min(1, max(histogram));
bgConfInitLB(4) = 0;

bgConfInitUB = [];
bgConfInitUB(1) = length(histogram);
bgConfInitUB(2) = max(length(histogram)/2, sigmaInit);
bgConfInitUB(3) = 2*max(histogram);
bgConfInitUB(4) = max(histogram)/2;

% fit gaussian
bgConf = lsqcurvefit(@fitGaussian1D, bgConfInit, length(histogram), histogram'...
	, bgConfInitLB, bgConfInitUB, optimset('Display', 'off'));

% as sigma is used always squared the algorithm might find it to be negativ. in this case we take
% the absolute value
bgConf(2) = abs(bgConf(2));

if debug == 1
	% display gauss distribution
	figure,plot(gaussian1D(bgConf, [1 1 1], length(histogram)));
end

% create binary mask for the Foreground of the image
FGMask = im;
FGMask(FGMask<bgConf(1)+bgWidth*bgConf(2)/2.3548200450309) = 0;
FGMask(FGMask>0) = 1;

% perform morphological opening to fill gaps and to delete noise
SE = ones(zeros(1, length(size(im)))+length(size(im)));
FGMask = imopen(FGMask, SE);
% FGMask = bwmorph(FGMask, 'fill');

if (debug == 1) && (length(size(im))==2)
	% display only Foreground
	figure,imshow(FGMask, []);
end

% return mean and standard deviation
bgMean = bgConf(1);
bgRMS = bgConf(2)/2.3548200450309;

end % function FGMask = segmentBackground3D(im, sigmaInterval, debug)

%% subfunction gauss1D = fitGaussian1D (gConf, sizex)
function gauss1D = fitGaussian1D(gConf, sizex)

gauss1D = gaussian1D(gConf, 1, sizex);

end % function gauss1D = gaussFit1D(gConf, sizex)
