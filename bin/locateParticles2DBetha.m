%
% function [gConfs, imMarked] = locateParticles2DBetha (im, dv, peakThreshRel, synapseDistInterval, gType, preprocess, displayResults, CCDSensitivity, EMGain)
%
% locates particles within one image
%
% parameters:
% -im			- the images containing the particles
% -dv			- real dimensions of one voxel 
% -tracking params	- set of parameters that might be necessary for tracking
%					NOTE: this parameter is a cell array!!
%					NOTE: each program picks the one it requires thus sometimes unused parametrs are
%					handed over
%		[01] - maxVoxelDisplacement -> 3
% 		[02] - synapseDistInterval -> 1
% 		[03] - minMaxFWHM -> [200, 200; 1400 1400]
% 		[04] - gType -> 1
% 		[05] - peakThreshRel -> 2
% 		[06] - minNumPhotonsPerPeak -> 25
% 		[07] - fitQualityThresholdFactor -> 2
% 		[08] - gaussianPart2Fit -> 2
%		[09] - intensityOutliers -> 0.1
%		[10] - reasonabilityChecks -> [1 1 1 1 1]
%		[11] - subImSizeFactor -> 2
%		[12] - useWatershedTransformRegions -> 1
%		[13] - weightedLsqPoisson -> [] (empty matrix means don't use weighted LSQ, put a number for adjusting, so far the number doesn't matter)
%		[14] - initFwhm -> [1200 1200] fwhm values in nm used to initialize the fit of the gaussian
%		[15] - watershedRegionControlAxisSizeFactor -> 1.5
%		[16] - gaussianOffsetConstraints -> 1
%		[17] - removeAllParticlesTooCloseToOthers -> 0
%		[18] - useImTopHatFilter -> [3 11 17]
%		[19] - doEMFitting -> 0 NOTE: if set to >0 than the useWatershedTransformRegions flag is ignored and treated as 0 thus a certain rectangular region is
%		used for fitting
%		[20] - useConfForSPT -> if set to 0 then all is as usual and the definition of parameters makes sense for detecting synapses with varying sizes
%								if set to 1 then the code is optimized for SPT where all particles are basically similar and where the definition of distances
%								and the size of the fitting regions are set constant and do not depend on the size of a particle anymore.
%								thus the parameters: 2 (synapseDistInterval) is now interpreted as the absolute distance in pixels that two particles must have
%								and 11 (subImSizeFactor) is an absolute value in pixels for the size of the fitting region
%		[21] - doSimpleButRobustThresholding -> if set to 1 than the new thresholding using the median is used and also the median is used for the
%		imTopHatFilter
%		[22] - useHighResolutionTopHatFilter -> if set to 1 than the new topHatFiler is used that is more robust (by using the mean instead of max) and allows to distinguish nearby synapses (NOTE: its slower)
%		trackingParams{1}=3,trackingParams{2}=1,trackingParams{3}=[200, 200; 1400 1400],trackingParams{4}=1,trackingPara
%		ms{5}=2,trackingParams{6}=25,trackingParams{7}=2,trackingParams{8}=2,trackingParams{9}=0.1
%		, trackingParams{10}=[1 1 1 1 1], trackingParams{11}=2, trackingParams{12}=1,
%		trackingParams{13}=[]
% -preprocess	- 1 if the image shall be preprocessed, 0 otherwise
% -display		- 1 if debug and result information shall be created, 0 otherwise
% -CCDSensitivity	- A/D converter discretization step size
% -EMGain			- EMGain used
% -imMask		-logical mask that defines a region where localizations are allowed. empty mask means no mask
% -initialLocations		- provides an initial vector with xy indexes for peak indications (primarily used in combination with magnus' fast locate)
function [gConfs, imMarked, residuals] = locateParticles2DBetha (im, dv, trackingParams, preprocess, display, CCDSensitivity, EMGain, imMask, initialLocations)

%% constructor

residuals = [];

% compute number of dimensions (just that i don't have to use the number all time)
dim = length(dv);
if display
	% save the original image
	imOrig = im;
end
% convert image to double
im = cast(im, 'double');

%%% get tracking params
% factor for the size of the synapse within which no other synapse is allowed to occure
% if useConfForSPT==1 then it is the absolute distance in pixels that two particles must have
synapseDistInterval = trackingParams{2};
% minimal and maximum values in nm for the FWHM to be accepted for a gaussian fit. example: [minx miny;maxx maxy]
minMaxFWHM = trackingParams{3};
% type of the gaussian:
% NOTE: in combination with single EM this option is ignored and for EMGM only option 0 and 2 are accepted
% -gType		- (0) -> gaussian with possibly varying width in x and y direction
%				- (1) -> gaussian with the same width in x and y direction
%				- (2) -> elliptical gaussian having varying x and y widths but can also be turned by
%				an angle theta
gType = trackingParams{4};
% relative parameter to estimate a threshold for the intensity of the gaussian peaks that are
% included in the gaussian fit procedure this values is used as a factor for the standard deviation of the mean image intensity. the
% produkt is added to the mean image intensity and is used as the threshold for peaks
% NOTE: this is related to thompson2002
peakThreshRel = trackingParams{5};
% minimal number of photons that is required to accept the peak as a real one
minNumPhotonsPerPeak = trackingParams{6};
% factor for the estimated noise at a pixel where the product is used as the threshold for deciding
% if a fitted gaussian is a sufficiently good fit to the peak
fitQualityThresholdFactor = trackingParams{7};
% is a factor for the standard deviation of the fitted gaussian. all values inbetween this range are
% used to control the quality of the gaussian fit
gaussianPart2Fit = trackingParams{8};
% upper percentige of the intensities in the image that are ignored when computing the threshold
% for the peaks. this is used to account for outliers in the image
intensityOutliers = trackingParams{9};
% defines which of the reasonability checks shall be performed
reasonabilityChecks = trackingParams{10};
% factor for the size of the subimage taken to fit a gaussian peak
% if useConfForSPT==1 then it is an absolute value in pixels for the size of the fitting region
subImSizeFactor = trackingParams{11};
% decides weather the regions in the image that the Gaussian is fitted to is created by a watershed
% transform (1) or just from a rectangular region arround the local maximas
useWatershedTransformRegion = trackingParams{12};
% if defined than a weighting parameter adjusted for poisson distributions is used for the lsq agorithm
weightedLsqPoisson = trackingParams{13};
% fwhm values in nm used to initialize the fit of the gaussian
initFwhm = trackingParams{14};
% factor for the mean minimal allowed fwhm of a peak, where the product is used to check that the
% smallest axis length of a watershed region has a certain size such that the fit is reasonable
watershedRegionControlAxisSizeFactor = trackingParams{15};
% defines how the offset for the gaussian is to be constraint during fitting
% 1 - the gaussian offset is limited by the estimate of the background +- a factor of its standard deviation
% 2 - the gaussian offset remains fixed to the estimate of the background
% 3 - the gaussian offset is nearly unbounded thus allows to fit also maybe only the top of a peak (a little comparable to rogers2006 BUT without gaussian weighting)
gaussianOffsetConstraints = trackingParams{16};
% if set to 1 then both particles not only the darker particle which is within synapseDistInterval (NOTE: this is a relative value) are removed such that we
% fit only particles that have no disturbances within their fitting window
% NOTE: if we use watershed transform than this parameter should be kept 0
removeAllParticlesTooCloseToOthers = trackingParams{17};
% if it is defined then tophat filtering is applied to find synaptic candidates instead of pure regional max finding
% setting: useImTopHatFilter = [iDiffRel, dTop, dBrim] == [relative intensity difference (as a factor for the std), diameter of top and brim]
useImTopHatFilter = trackingParams{18};
% if set to >0 than EM fitting instead of least square curve fitting is used and the number defines the number of iterations to be used by EM
% NOTE: if set to >0 than the useWatershedTransformRegions flag is ignored and treated as 0 thus a certain rectangular region is used for fitting
doEMFitting = trackingParams{19};
% if set to 0 then all is as usual and the definition of parameters makes sense for detecting synapses with varying sizes if set to 1 then the code is
% optimized for SPT where all particles are basically similar and where the definition of distances and the size of the fitting regions are set constant
% and do not depend on the size of a particle anymore.
% thus the parameters: 2 (synapseDistInterval) is now interpreted as the absolute distance in pixels that two particles must have
% and 11 (subImSizeFactor) is an absolute value in pixels for the size of the fitting region
useConfForSPT = trackingParams{20};
% if set to 1 than the new thresholding using the median is used and also the median is used for the imTopHatFilter
doSimpleButRobustThresholding = trackingParams{21};
% if set to 1 than the new topHatFiler is used that is more robust (by using the mean instead of max) and allows to distinguish nearby synapses (NOTE: its slower)
useHighResolutionTopHatFilter = trackingParams{22};

%% parameters/variables

% standard deviation for the gaussian used for smoothing via convolution (crocker1996) and parameters from (rogers2007)
stdGaussianSmoothing = 1;

% threshold for the acceptance of a normal distribution of the residuals of a fit
fitAlphaThresh = 0.05;

% size of the partial image that is taken to fit the gauss
% NOTE: the values have to be odd!
% TODO: these values have to be determined automatically
% imPeakSize = [min(13, size(im, 2)), min(13, size(im, 1))];
if ~useConfForSPT
	imPeakSize = [min(round(minMaxFWHM(2, 1)*subImSizeFactor/dv(1)), size(im, 2)), min(round(minMaxFWHM(2, 2)*subImSizeFactor/dv(2)), size(im, 1))];
else
	% if localization is used for SPT than take a constant region
	imPeakSize = [subImSizeFactor, subImSizeFactor];
end
% correct imPeaksize if necessary
if mod(imPeakSize(1), 2) == 0
	imPeakSize(1)=imPeakSize(1)-1;
end
if mod(imPeakSize(2), 2) == 0
	imPeakSize(2)=imPeakSize(2)-1;
end

% half the size of the partial image that is taken to fit the gauss
imPeakSizeHalf = [floor(imPeakSize(1)/2), floor(imPeakSize(2)/2)];
% radius of the circle that acts as a marker for the peaks that were found
markerRadius = 5;
% interval in number of sigmas for the background noise
bgNoiseInterval = 3;
% interval in number of sigmas for the gaussian that is used to collect the number of photons of the
% peak
photonNumInterval = 3;

% preprocess the image if required
if preprocess == 1
	%%% smooth the data such that we can easily find the local maxima of the peak
	% h = ones(3, 3) / 9;
	% convolution of the image with a gaussian (crocker1996) and parameters from (rogers2007)
	h = fspecial('gaussian', [3 3], stdGaussianSmoothing);
	im = imfilter(im, h, 'replicate', 'same', 'conv');
end

% do some threshold computation according to thompson2002
if ~doSimpleButRobustThresholding
	% compute threshold for peak intensity (related to thompson2002)
	% peakThresh = 0.65*maxPeak
	% peakThresh = peakThreshRel*std(im(:))+mean(im(:))
	% substract the background from the image intensities
	[null bgNoiseMean bgNoiseRMS] = segmentBackground(round(im), bgNoiseInterval, 0);

	signalValues = im-(bgNoiseMean+bgNoiseInterval*bgNoiseRMS);
	% get only values above background
	signalValues = signalValues(signalValues>0);

	% ignore the upper x percent of the image values to avoid influence of outliers
	% thus sort
	signalValues = sort(signalValues);
	% and then remove the upper part
	signalValues = signalValues(1:length(signalValues)-floor(length(signalValues)*intensityOutliers));

	% compute threshold
	peakThresh = peakThreshRel*std(signalValues)+mean(signalValues)+(bgNoiseMean+bgNoiseInterval*bgNoiseRMS);
% do simple but robust thresholding using median
else
	bgNoiseMean = median(im(:));
	bgNoiseRMS = mad(im(:), 1);
	peakThresh = bgNoiseMean+peakThreshRel*bgNoiseRMS;
end


% matrix containing the configurations of all gaussian peaks that have been detected in the image
gConfs = [];

%%% default initial configuration for the gaussian fit
% TODO: either classify the configuration or let it be derived from the image
gConfInit = [];
% x position of the center
gConfInit(1) = ceil(imPeakSize(1)/2);
% y position of the center
gConfInit(2) = ceil(imPeakSize(2)/2);
% full width x at half maximum (fwhm) of the gauss
gConfInit(3) = initFwhm(1)/dv(1);%(minMaxFWHM(2, 1)/2)/dv(1); % 1.7;
% full width y at half maximum (fwhm) of the gauss
gConfInit(4) = initFwhm(2)/dv(2);%(minMaxFWHM(2, 2)/2)/dv(2); % 1.7;
% indicator for the height of the gaussian peak
gConfInit(5) = peakThresh-bgNoiseMean; % 1;
% offset for the gauss
gConfInit(6) = bgNoiseMean; % 0;
% % real x position of the center in the image
% gConfInit(7) = -1;
% % real y position of the center in the image
% gConfInit(8) = -1;
% if we use most general gaussian than we need the angle theta as well
if gType == 2
	gConfInit(9) = 0;
end

%%% lower bound for the configuration of the gaussians that are fitted to the peaks
% TODO: either classify the configuration or let it be derived from the image
gConfLowerBound = [];
% x position of the center
gConfLowerBound(1) = gConfInit(1)-imPeakSizeHalf(1)+1;
% y position of the center
gConfLowerBound(2) = gConfInit(2)-imPeakSizeHalf(2)+1;
% full width x at half maximum (fwhm) of the gauss
% NOTE: the minimum for the fit is smaller than the allowed size to be able to discard it if it would
% be too small
gConfLowerBound(3) = (minMaxFWHM(1, 1)/2)/dv(1); % 0.0001;
% full width y at half maximum (fwhm) of the gauss
% NOTE: the minimum for the fit is smaller than the allowed size to be able to discard it if it would
% be too small
gConfLowerBound(4) = (minMaxFWHM(1, 2)/2)/dv(2); % 0.0001;
% indicator for the height of the gaussian peak
% TODO: this should be computed from the intensity threshold for the peaks
gConfLowerBound(5) = 1;
% offset for the gauss
if gaussianOffsetConstraints == 1
	gConfLowerBound(6) = max(1, bgNoiseMean-bgNoiseInterval*bgNoiseRMS);
elseif gaussianOffsetConstraints == 2
	gConfLowerBound(6) = gConfInit(6)-0.0001;
elseif gaussianOffsetConstraints == 3
	gConfLowerBound(6) = 0;
end
% if we use most general gaussian than we need the angle theta as well
if gType == 2
	gConfLowerBound(9) = -10*360;
end


%%% upper bound for the configuration of the gaussians that are fitted to the peaks
% TODO: either classify the configuration or let it be derived from the image
gConfUpperBound = [];
% x position of the center
gConfUpperBound(1) = gConfInit(1)+imPeakSizeHalf(1)-1;
% y position of the center
gConfUpperBound(2) = gConfInit(2)+imPeakSizeHalf(2)-1;
% full width x at half maximum (fwhm) of the gauss
% NOTE: the maximum for the fit is higher than the allowed size to be able to discard it if it would
% be too big
gConfUpperBound(3) = (minMaxFWHM(2, 1)*1.5)/dv(1); % max(imPeakSize);
% full width y at half maximum (fwhm) of the gauss
% NOTE: the maximum for the fit is higher than the allowed size to be able to discard it if it would
% be too big
gConfUpperBound(4) = (minMaxFWHM(2, 2)*1.5)/dv(2); % max(imPeakSize);
% indicator for the height of the gaussian peak
% TODO: check if this makes sense
gConfUpperBound(5) = 1.2*(max(max(im))-bgNoiseMean);%2*max(max(im));
% offset for the gauss
% TODO: maybe a dependency on the SNR makes sense
if gaussianOffsetConstraints == 1
	% % NOTE: isn't that a mistake that i use max here? shouldn't it be min!!
	gConfUpperBound(6) = bgNoiseMean+bgNoiseInterval*bgNoiseRMS;%max(max(max(im))/2, bgNoiseMean+bgNoiseInterval*bgNoiseRMS);
elseif gaussianOffsetConstraints == 2
	gConfUpperBound(6) = gConfInit(6)+0.0001;
elseif gaussianOffsetConstraints == 3
	gConfUpperBound(6) = max(max(im));
end

% if we use most general gaussian than we need the angle theta as well
if gType == 2
	gConfUpperBound(9) = 10*360;
end

%% preliminaries

% load EM library
if doEMFitting
	savedWarningState = warning ('off','all');
% 	loadlibrary('/home/stefan/PointillismEM/mdaEM/pem.so', '/home/stefan/PointillismEM/mdaEM/pem.h');
% 	loadlibrary('/home/stefan/PointillismEM_twoGaussian/mdaEM/pem2.so', '/home/stefan/PointillismEM_twoGaussian/mdaEM/pem2.h');
	warning(savedWarningState);
end

%%% try to make two very close synapses distinguishable
% usually it happens that two close synapses are fitted as one together. this shall be avoided via
% creating areas arround the local maxima (using watershed) such that each synapse is taken as its
% own
if useWatershedTransformRegion && ~doEMFitting
	% invert the image first so that minima become maxima and vice versa
% 	imWS = cast(abs(im-(max(max(im)))), 'uint16');
	imWS = cast(max(max(im))-im, 'uint16');
	% do watersheed
	imWS = watershed(imWS, 8);

% 	% plot watershed result
% 	figure, imshow(label2rgb(imWS));
end

%%% detect possible peaks that will be fitted later on

% if no initial locations have been provided than find some in the image
if isempty(initialLocations)

	% find local maximas in the image and save as binary
	if isempty(useImTopHatFilter)
		% find all local maxima
		imLocMax = imregionalmax(im, 8);
	else
	% 	stdSignalValues = std(signalValues)
	% 	peakThresh
		% find local maxima using top hat filtering
		if ~useHighResolutionTopHatFilter
			if ~doSimpleButRobustThresholding
				imLocMax = topHatFilter(im, std(signalValues).*useImTopHatFilter(1), useImTopHatFilter(2), useImTopHatFilter(3), 1 , 1);
			else
				imLocMax = topHatFilter(im, bgNoiseRMS.*useImTopHatFilter(1), useImTopHatFilter(2), useImTopHatFilter(3), 1 , 1);
			end
		else
			if ~doSimpleButRobustThresholding
				imLocMax = topHatFilterHighResolution(im, std(signalValues).*useImTopHatFilter(1), useImTopHatFilter(2), useImTopHatFilter(3), 1);
			else
				imLocMax = topHatFilterHighResolution(im, bgNoiseRMS.*useImTopHatFilter(1), useImTopHatFilter(2), useImTopHatFilter(3), 1);
			end
		end
	end

	% get positions of all local maxima
	% [locMaxY, locMaxX] = find(imLocMax==1);
	locMaxX = [];
	locMaxY = [];
	for x=1:size(imLocMax, 2)
		for y=1:size(imLocMax, 1)
			if imLocMax(y, x) == 1
				locMaxX = [locMaxX; x];
				locMaxY = [locMaxY; y];
			end
		end
	end
else
	% create location index from initial positions
	locMaxX = initialLocations(:, 1);
	locMaxY = initialLocations(:, 2);
	
	% if we do EM fitting using gaussian mixtures we create a map of the local maxima because if we get the peaks from outside it doesn't exist, yet
	if doEMFitting
		imLocMax = uint16(zeros(size(im)));
		for i=1:length(locMaxX)
			imLocMax(locMaxY(i), locMaxX(i)) = 1;
		end
	end
end

%% plot the results of tophatfilter and relative thresholding
% figure, imagesc(im), colormap(gray), hold on, plot(locMaxX, locMaxY, 'b.', 'LineStyle', 'none')
% locMaxXTmp = locMaxX;
% locMaxYTmp = locMaxY;
% for i = fliplr(1:length(locMaxXTmp))
% 	if im(locMaxYTmp(i), locMaxXTmp(i)) <= peakThresh
% 		locMaxXTmp(i) = [];
% 		locMaxYTmp(i) = [];
% 	end
% end
% plot(locMaxXTmp, locMaxYTmp, 'r.', 'LineStyle', 'none'), hold off;

%% main
%% detect gaussian peaks

% init number of peaks
peakNum = 0;

% loop over all local maximas
for i=1:length(locMaxX)
	
	% only accept a local maxima as a peak if it is above the peak threshold
	if im(locMaxY(i), locMaxX(i)) > peakThresh && (isempty(imMask) || imMask(locMaxY(i), locMaxX(i))) %%&& (locMaxX(i)>10) && (locMaxY(i)>10)

		% increase number of peaks
		peakNum = peakNum+1;

		% get the initial configuration for the peak
		gConfInitPeak = gConfInit;

		% cut the corresponding part in the image either using a watershed transform region or taking a
		% rectangular region arround the maxima
		if ~useWatershedTransformRegion || (useWatershedTransformRegion && doEMFitting)
			% use rectangular region
			
			% compute the indexes of the peakImage in the real image
			imPeakIdx = zeros(2, 2);

			imPeakIdx(1, 1) = locMaxX(i)-imPeakSizeHalf(1);
			imPeakIdx(1, 2) = locMaxX(i)+imPeakSizeHalf(1);
			imPeakIdx(2, 1) = locMaxY(i)-imPeakSizeHalf(2);
			imPeakIdx(2, 2) = locMaxY(i)+imPeakSizeHalf(2);

			% correct the peak index and the initial configuration if the peak is so close to any border
			% that not the full peakImage with the peak in the center can be obtained
			% NOTE: in every dimension the border can be crossed at most in one direction as the size of
			% the peakImage is known to be at most the size of the real image

			% if over left border (x)
			if imPeakIdx(1, 1)<1
				imPeakIdx(1, 2) = imPeakIdx(1, 2)+(1+abs(imPeakIdx(1, 1)));
				gConfInitPeak(1) = gConfInitPeak(1)-(1+abs(imPeakIdx(1, 1)));
				imPeakIdx(1, 1) = 1;

			% if over right border (x)
			elseif imPeakIdx(1, 2)>size(im, 2)
				imPeakIdx(1, 1) = imPeakIdx(1, 1)-(imPeakIdx(1, 2)-size(im, 2));
				gConfInitPeak(1) = gConfInitPeak(1)+(imPeakIdx(1, 2)-size(im, 2));
				imPeakIdx(1, 2) = size(im, 2);
			end
			% if over upper border (y)
			if imPeakIdx(2, 1)<1
				imPeakIdx(2, 2) = imPeakIdx(2, 2)+(1+abs(imPeakIdx(2, 1)));
				gConfInitPeak(2) = gConfInitPeak(2)-(1+abs(imPeakIdx(2, 1)));
				imPeakIdx(2, 1) = 1;
			% if over lower border (y)
			elseif imPeakIdx(2, 2)>size(im, 1)
				imPeakIdx(2, 1) = imPeakIdx(2, 1)-(imPeakIdx(2, 2)-size(im, 1));
				gConfInitPeak(2) = gConfInitPeak(2)+(imPeakIdx(2, 2)-size(im, 1));
				imPeakIdx(2, 2) = size(im, 1);
			end

			% finally get the corresponding part in the image, which is our peak where a gauss shall be fitted to
			imPeak = im(imPeakIdx(2,1):imPeakIdx(2, 2)...
				, imPeakIdx(1,1):imPeakIdx(1, 2));

			null = gConfUpperBound-gConfLowerBound;
			null(null~=0) = 1;
			if sum(null, 2)<6
				gConfUpperBound
				gConfLowerBound
			end

			% if we want to do EM fitting instead of LSQ fitting
			if doEMFitting

				%% EM fitting with gaussian mixture of multiple elliptical gaussian

				% check if the fit shall be symmetric or not
				if gType~=0 && gType~=2
					fprintf('ERROR: incorrect gType option (%d) only 1 and 2 are accepted options - exiting\n', gType);
					return;
				end
				
				%%% find all peak indications in the image such that EM can be initialized
				
				% also cut the region from the image containing the local maxima
				imPeakLocMax = imLocMax(imPeakIdx(2,1):imPeakIdx(2, 2)...
				, imPeakIdx(1,1):imPeakIdx(1, 2));

				% test that the current peak indication exists in the peak image
				if imPeakLocMax(gConfInitPeak(2), gConfInitPeak(1)) ~= 1
					fprintf('ERROR: peak indication does not exist in peak image - exiting\n');
					return;
				end
				
				% remove current peak from local peak image such that we do not find it again
				imPeakLocMax(gConfInitPeak(2), gConfInitPeak(1)) = 0;
				
% 				figure, imshow(imPeakLocMax,  []);
				
				% get positions of current and all other local maxima (!!!!!!!!that pass the threshold!!!!!!!)
				peakLocMax = [gConfInitPeak(1) gConfInitPeak(2)];
				for x=1:size(imPeakLocMax, 2)
					for y=1:size(imPeakLocMax, 1)
						if imPeakLocMax(y, x) == 1 && imPeak(y, x) > peakThresh
							peakLocMax = [peakLocMax; x y];
						end
					end
				end
				
% 				%%% loop over possible number of peaks in the region and fit
% 				for numPeaksInit = size(peakLocMax, 1)
				% NOTE: so far we do not loop but simply trust on the initialization
				numPeaksInit = size(peakLocMax, 1);
				
					% create initialization
					% [number of background photons, number of photons in each gaussian (thus n values), for each gaussian: MuX, MuY, SigmaX, SigmaY, (Rho)]
% 					% NOOOO NOTE: for initialization x and y are provided in the order of matlabs matrix definition BUT the output is again in the opposite definition, thus, x before y
					
% % 					gmConfInit = 1/(numPeaksInit+1);
% % 					gmConfInit = ones(1, numPeaksInit+1).*gmConfInit;
					gmConfInit = [gConfInitPeak(6)*numel(imPeak) ones(1, numPeaksInit).*(gConfInitPeak(5)*10)];
% 					gmConfInit = ones(1, numPeaksInit+1);
					for numPeaksInitIdx = 1:numPeaksInit
% 						gmConfInit = [gmConfInit ceil(imPeakSize(2)/2) ceil(imPeakSize(1)/2) initFwhm(2)/dv(2)/2.3548200450309 initFwhm(1)/dv(1)/2.3548200450309 0];

						% if symmetrical Gaussians shall be fitted
						if gType == 0
							gmConfInit = [gmConfInit peakLocMax(numPeaksInitIdx, 1:2) initFwhm(1)/dv(1)/2.3548200450309 initFwhm(2)/dv(2)/2.3548200450309];
% 							gmConfInit = [gmConfInit peakLocMax(numPeaksInitIdx, 2) peakLocMax(numPeaksInitIdx, 1) initFwhm(1)/dv(1)/2.3548200450309 initFwhm(2)/dv(2)/2.3548200450309];
						% else if also the angle rho can be adjusted
						else
							gmConfInit = [gmConfInit peakLocMax(numPeaksInitIdx, 1:2) initFwhm(1)/dv(1)/2.3548200450309 initFwhm(2)/dv(2)/2.3548200450309 0];
% 							gmConfInit = [gmConfInit peakLocMax(numPeaksInitIdx, 2) peakLocMax(numPeaksInitIdx, 1) initFwhm(1)/dv(1)/2.3548200450309 initFwhm(2)/dv(2)/2.3548200450309 0];
						end
					end
					
					% do EM fitting
					% if symmetrical Gaussians shall be fitted
					if gType == 0
						[gmConf, imModel] = PointillismEMmex(imPeak, doEMFitting, gmConfInit, 1);
					% else if also the angle rho can be adjusted
					else
						[gmConf, imModel] = PointillismEMmex(imPeak, doEMFitting, gmConfInit, 0);
					end
% 					gmConf
					% loop over gConf and adjust for the offset in the coordinate system of EMGM
					% NOTE: the difference stems from the fact that EMGM really integrates over the pixels so e.g. all photons are at pixel on that are >0 and
					% <=1. however when computing a surface pixel 1 means position 1 but in EMGM pixel 1 actually means value 0.5 because at this position +-0.5
					% all photons were integrated
					for numPeaksInitIdx = 1:numPeaksInit
						% if symmetrical Gaussians shall be fitted
						if gType == 0
							gmConf(1+numPeaksInit+(numPeaksInitIdx-1)*4+1) = gmConf(1+numPeaksInit+(numPeaksInitIdx-1)*4+1)+0.5;
							gmConf(1+numPeaksInit+(numPeaksInitIdx-1)*4+2) = gmConf(1+numPeaksInit+(numPeaksInitIdx-1)*4+2)+0.5;
						% else if also the angle rho can be adjusted
						else
							gmConf(1+numPeaksInit+(numPeaksInitIdx-1)*5+1) = gmConf(1+numPeaksInit+(numPeaksInitIdx-1)*5+1)+0.5;
							gmConf(1+numPeaksInit+(numPeaksInitIdx-1)*5+2) = gmConf(1+numPeaksInit+(numPeaksInitIdx-1)*5+2)+0.5;
						end
					end

					% gmConf definition:
					% [number of background photons, number of photons in each gaussian (thus n values), MuX, MuY, SigmaX, SigmaY, (rho)]

					% create something like gConfs BUT extract only the information of the first peak
%					-(1)	- center x position
%					-(2)	- center y position
%					-(3)	- fwhm x, gaussian with major axis in x direction
%					-(4)	- fwhm y, gaussian with major axis in x direction
%				    -(5)	- height of the gaussian with major axis in x direction
%				    -(6)	- offset
%				    .. reserved
%					-(9) - turning angle theta if an elliptical gaussian is used
%					-..	- further parameters are up to you but are ignored here
					% if symmetrical Gaussians shall be fitted
					if gType == 0
						null = [gmConf(1+numPeaksInit+1) gmConf(1+numPeaksInit+2) gmConf(1+numPeaksInit+3)*2.3548200450309 gmConf(1+numPeaksInit+4)*2.3548200450309...
							gmConf(2) gmConf(1) gmConf(1+numPeaksInit+1) gmConf(1+numPeaksInit+2)];
					% else if also the angle rho can be adjusted
					else
						null = [gmConf(1+numPeaksInit+1) gmConf(1+numPeaksInit+2) gmConf(1+numPeaksInit+3)*2.3548200450309 gmConf(1+numPeaksInit+4)*2.3548200450309...
							gmConf(2) gmConf(1) gmConf(1+numPeaksInit+1) gmConf(1+numPeaksInit+2) gmConf(1+numPeaksInit+5)];
					end
% 					null
					
					%% draw fitting results
% 					figure, imagesc(imPeak), colormap(gray), hold on;
% 					% draw initializations as points
% 					plot(peakLocMax(:, 1), peakLocMax(:, 2), 'g*');
% 					
% 					if ~isnan(gmConf(1, 1))
% 						for numPeaksInitIdx = 1:numPeaksInit
% 							% if symmetrical Gaussians shall be fitted
% 							if gType == 0
% 								if numPeaksInitIdx == 1
% 									drawGaussian([null(1) null(2) null(3)/2.3548200450309 null(4)/2.3548200450309 0], 'r');
% 								else
% 									drawGaussian([gmConf(1+numPeaksInit+(numPeaksInitIdx-1)*4+1) gmConf(1+numPeaksInit+(numPeaksInitIdx-1)*4+2)...
% 										gmConf(1+numPeaksInit+(numPeaksInitIdx-1)*4+3) gmConf(1+numPeaksInit+(numPeaksInitIdx-1)*4+4)...
% 										0], 'b');
% 								end
% 							% else if also the angle rho can be adjusted
% 							else
% 								if numPeaksInitIdx == 1
% 									drawGaussian([null(1) null(2) null(3)/2.3548200450309 null(4)/2.3548200450309 null(9)], 'r');
% 								else
% 									drawGaussian([gmConf(1+numPeaksInit+(numPeaksInitIdx-1)*5+1) gmConf(1+numPeaksInit+(numPeaksInitIdx-1)*5+2)...
% 										gmConf(1+numPeaksInit+(numPeaksInitIdx-1)*5+3) gmConf(1+numPeaksInit+(numPeaksInitIdx-1)*5+4)...
% 										gmConf(1+numPeaksInit+(numPeaksInitIdx-1)*5+5)], 'b');
% 								end
% 							end
% 						end
% 					end
% 					hold off;

					% if the fitting didn't yield a result than skip
					if isnan(null(1))
						% decrease number of found particles
						peakNum = peakNum-1;
						% jump to the end of the loop
						continue;
					end
% 				end

				%% EM fitting with single elliptical gaussian
				
% 				% create initialization
% 				% MuX, MuY, SigmaX, SigmaY, Rho, Background fraction
% 				gmConfInit = [ceil(imPeakSize(1)/2) ceil(imPeakSize(2)/2) initFwhm(1)/dv(1)/2.3548200450309 initFwhm(2)/dv(2)/2.3548200450309 0 0.5];
% 
% 				gmConf = calllib('pem', 'pem', imPeak, gmConfInit, doEMFitting, 1e-10)';
% 
% 				% create something like gConfs
% % 			    -(1)	- center x position
% % 				-(2)	- center y position
% % 				-(3)	- fwhm x, gaussian with major axis in x direction
% % 				-(4)	- fwhm y, gaussian with major axis in x direction
% %  			    -(5)	- height of the gaussian with major axis in x direction
% % 			    -(6)	- offset
% % 			    .. reserved
% % 				-(9) - turning angle theta if an elliptical gaussian is used
% % 				-..	- further parameters are up to you but are ignored here
% 				null = [gmConf(1:2) gmConf(3)*2.3548200450309 gmConf(4)*2.3548200450309 1-gmConf(6) gmConf(6) gmConf(1:2) gmConf(5)];
% 
% 				% draw fitting results
% 				figure, imagesc(imPeak), colormap(gray), hold on;
% 				% draw initializations as points
% 				plot(gmConfInit(1), gmConfInit(2), 'g*');
% 				if ~isnan(null(1, 1))
% 					drawGaussian([null(1) null(2) null(3)/2.3548200450309 null(4)/2.3548200450309 null(9)], 'r');
% 				end
% 				hold off;

				%% EM fitting with sum of two elliptical gaussian
				
% 				% create initialization
% 				% MuX, MuY, MuX2, MuY2, SigmaX, SigmaY, SigmaX2, SigmaY2, Rho, Rho2, Background fraction
% 				gmConfInit = [ceil(imPeakSize(1)/2)+1 ceil(imPeakSize(2)/2)+1 ceil(imPeakSize(1)/2)-1 ceil(imPeakSize(2)/2)-1 ...
% 					initFwhm(1)/dv(1)/2.3548200450309 initFwhm(2)/dv(2)/2.3548200450309 initFwhm(1)/dv(1)/2.3548200450309 initFwhm(2)/dv(2)/2.3548200450309...
% 					0 0 0.5];
% 
% 				gmConf = calllib('pem2', 'pem2', imPeak, gmConfInit, doEMFitting, 1e-10)';
% 
% 				% create something like gConfs
% % 			    -(1)	- center x position
% % 				-(2)	- center y position
% % 				-(3)	- fwhm x, gaussian with major axis in x direction
% % 				-(4)	- fwhm y, gaussian with major axis in x direction
% %   			-(5)	- height of the gaussian with major axis in x direction
% % 			    -(6)	- offset
% % 			    .. reserved
% % 				-(9) - turning angle theta/rho if an elliptical gaussian is used
% % 			    .. reserved
% % 				-(12) 	- height of gaussian with major axis in y direction
% % 			    -(13)	- fwhm x, gaussian with major axis in y direction
% % 			    -(14)	- fwhm y, gaussian with major axis in y direction
% % 			    -(15)	- center x position, gaussian with major axis in y direction
% % 			    -(16)	- center y position, gaussian with major axis in y direction
% % 				-(17) - turning angle theta/rho2
% % 				-..	- further parameters are up to you but are ignored here
% 				null = [gmConf(1:2) gmConf(5)*2.3548200450309 gmConf(6)*2.3548200450309 1-gmConf(11) gmConf(11) gmConf(1:2) gmConf(9)...
% 					nan nan  1-gmConf(11) gmConf(7)*2.3548200450309 gmConf(8)*2.3548200450309 gmConf(3:4) gmConf(10)];

			% else do normal LSQ curve fitting
			else
				% fit a gaussian to the peak in the image
				[null residuals] = fitGaussian2D(imPeak, dv, gConfInitPeak, gConfLowerBound, gConfUpperBound, gType, weightedLsqPoisson);

% 				% fit a real gaussian mixture to the peak in the image
% 				gConfInit(5) = max(1.1, gConfInit(5)/2);
% 				gConfLowerBound(5) = 1;
% 				gConfLowerBound(6) = 1;
% 				
% 				gConfInit(12) = gConfInit(5);
% 				gConfLowerBound(12) = gConfLowerBound(5);
% 				gConfUpperBound(12) = gConfUpperBound(5);
% 				gConfInit(13:14) = [gConfInit(4) gConfInit(3)];
% 				gConfLowerBound(13:14) = [gConfLowerBound(4) gConfLowerBound(3)];
% 				gConfUpperBound(13:14) = [gConfUpperBound(4) gConfUpperBound(3)];
% 				[null residuals] = fitGaussianMixture2D(imPeak, dv, gConfInit, gConfLowerBound, gConfUpperBound, weightedLsqPoisson);

% 				% fit our half gaussian mixture to the peak in the image
% 				gConfInit(5) = max(1.1, gConfInit(5)/2);
% 				gConfLowerBound(5) = 1;
% 				gConfLowerBound(6) = 1;
% 				
% 				gConfInit(12) = gConfInit(5);
% 				gConfInit(5) = max(1.1, gConfInit(5).*0.75);
% 				gConfLowerBound(12) = gConfLowerBound(5);
% 				gConfUpperBound(12) = gConfUpperBound(5);
% 				[null residuals] = fitGaussianMixture2DHalfDiscrete(im, dv, gConfInit, gConfLowerBound, gConfUpperBound, weightedLsqPoisson);

% 				% fit our gaussian sum where only the intensity between the two gaussian varies to the peak in the image
% 				gConfInit(5) = max(1.1, gConfInit(5)/2);
% 				gConfLowerBound(5) = 1;
% 				gConfLowerBound(6) = 1;
% 				
% 				gConfInit(12) = gConfInit(5);
% 				gConfInit(5) = max(1.1, gConfInit(5).*0.75);
% 				gConfLowerBound(12) = gConfLowerBound(5);
% 				gConfUpperBound(12) = gConfUpperBound(5);
% 				[null residuals] = fitGaussianSum2D(imPeak, dv, gConfInit, gConfLowerBound, gConfUpperBound, weightedLsqPoisson);

% 				% draw fitting results
% 				figure, imagesc(imPeak), colormap(gray), hold on;
% 				% draw initializations as points
% 				plot(gConfInit(1), gConfInit(2), 'g*');
% 				if ~isnan(null(1, 1))
% 					if gType == 2
% 						drawGaussian([null(1) null(2) null(3)/2.3548200450309 null(4)/2.3548200450309 null(9)], 'r');
% 					else
% 						drawGaussian([null(1) null(2) null(3)/2.3548200450309 null(4)/2.3548200450309 0], 'r');
% 					end
% 				end
% 				hold off;

			end

			% update the real xy positions of the gauss center in the image
			null(7) = null(1) -1 + imPeakIdx(1, 1);
			null(8) = null(2) -1 + imPeakIdx(2, 1);

			gConfs = [gConfs; null];
% 			gConfs = [gConfs; fitGaussian2D(imPeak, dv, gConfInitPeak, gConfLowerBound, gConfUpperBound, gType, weightedLsqPoisson)];

% 			% update the real xy positions of the gauss center in the image
% 			gConfs(peakNum, 7) = gConfs(peakNum, 1) -1 + imPeakIdx(1, 1);
% 			gConfs(peakNum, 8) = gConfs(peakNum, 2) -1 + imPeakIdx(2, 1);
		else
			% use watershed transform defined region
			
% 			% so far just warn if we hit only the border of a region
% 			if imWS(locMaxY(i), locMaxX(i))==0
% 				fprintf('WARNING: watershed transform border was found as a local maxima, so far unhandled!\n');
% 			end

			% if we hit only the border of a region with a local maxima than we ignore this maxima
			if imWS(locMaxY(i), locMaxX(i))==0
				fprintf('WARNING: watershed transform border was found as a local maxima - deleting this maxima!\n');

				% decrease number of found particles
				peakNum = peakNum-1;
				% jump to the end of the loop
				continue;
			end

			% adjust the position parameters of the configuration
			gConfInitPeak(1) = locMaxX(i);
			gConfInitPeak(2) = locMaxY(i);
			% adjust also the boundaries
			gConfLowerBoundPeak = gConfLowerBound;
			gConfUpperBoundPeak = gConfUpperBound;
			gConfLowerBoundPeak(1) = max(1, gConfInitPeak(1)-imPeakSizeHalf(1)+1);
			gConfLowerBoundPeak(2) = max(1, gConfInitPeak(2)-imPeakSizeHalf(2)+1);
			gConfUpperBoundPeak(1) = min(size(im, 2), gConfInitPeak(1)+imPeakSizeHalf(1)-1);
			gConfUpperBoundPeak(2) = min(size(im, 1), gConfInitPeak(2)+imPeakSizeHalf(2)-1);

			% take the full image but with only the according area of the watershed region remaining
			% as the real image information. all other values are set to nan
			imPeak = im;
			imPeak = imPeak.*(imWS==imWS(locMaxY(i), locMaxX(i)));
			imPeak(imPeak==0) = NaN;

			null = gConfUpperBoundPeak-gConfLowerBoundPeak;
			null(null~=0) = 1;
			if sum(null, 2)<6
				gConfUpperBoundPeak
				gConfLowerBoundPeak
			end
			
			% check if the number of pixels that the peak consists of is high enough for an lsq
			% fitting to work
			if sum(sum(~isnan(imPeak))) < length(gConfInitPeak)
				% decrease number of found particles
				peakNum = peakNum-1;
				% jump to the end of the loop
				continue;
			end
			
			if gConfLowerBoundPeak(1) >= gConfUpperBoundPeak(1) || gConfLowerBoundPeak(2) >= gConfUpperBoundPeak(2) || gConfLowerBoundPeak(3) >= gConfUpperBoundPeak(3)...
					|| gConfLowerBoundPeak(4) >= gConfUpperBoundPeak(4) || gConfLowerBoundPeak(5) >= gConfUpperBoundPeak(5) || gConfLowerBoundPeak(6) >= gConfUpperBoundPeak(6)
				
				fprintf('WARNING: computed upper and lower bounds are inconsumerate - deleting this maxima!\n');

				gConfInitPeak
				gConfLowerBoundPeak
				gConfUpperBoundPeak

				% decrease number of found particles
				peakNum = peakNum-1;
				% jump to the end of the loop
				continue;
			end
			
			% fit a gaussian to the peak in the image
			[null residuals] = fitGaussian2DDiscrete(imPeak, dv, gConfInitPeak, gConfLowerBoundPeak, gConfUpperBoundPeak, gType, weightedLsqPoisson);

% 			% fit a gaussian mixture to the peak in the image
% 			gConfInitPeak(5) = max(1.1, gConfInitPeak(5)/2);
% 			gConfLowerBoundPeak(5) = 1;
% 			gConfLowerBoundPeak(6) = 1;
% 			
% 			gConfInitPeak(12) = gConfInitPeak(5);
% 			gConfLowerBoundPeak(12) = gConfLowerBoundPeak(5);
% 			gConfUpperBoundPeak(12) = gConfUpperBoundPeak(5);
% 			[null residuals] = fitGaussianMixture2DDiscrete(im, dv, gConfInitPeak, gConfLowerBoundPeak, gConfUpperBoundPeak, weightedLsqPoisson);

			% update the real xy positions of the gauss center in the image
			% NOTE: because we take the full image (with just the watershed region left) coordinates
			% in 1,2 are the same as in 7,8
			null(7) = null(1);
			null(8) = null(2);

			gConfs = [gConfs; null];
% 			gConfs = [gConfs; fitGaussian2DDiscrete(imPeak, dv, gConfInitPeak, gConfLowerBoundPeak, gConfUpperBoundPeak, gType, weightedLsqPoisson)];
% 			gConfs = [gConfs; fitGaussian2DDiscreteWeighted(imPeak, dv, gConfInitPeak, gConfLowerBoundPeak, gConfUpperBoundPeak, gType, 1./cast(imPeak, 'double'))];
			
% 			% update the real xy positions of the gauss center in the image
% 			% NOTE: because we take the full image (with just the watershed region left) coordinates
% 			% in 1,2 are the same as in 7,8
% 			gConfs(peakNum, 7) = gConfs(peakNum, 1);
% 			gConfs(peakNum, 8) = gConfs(peakNum, 2);
			
			% TODO: sometimes the error: Warning: Cannot solve problems with fewer equations than variables and with bounds.
			% An error will be issued for this case in a future release. Ignoring bounds, using
			% Levenberg-Marquardt method.
			% may appear -> this has to be checked because unreasonable results may be the case
			% so far we just check the correct xy setting
			if round(gConfs(peakNum, 8)) > size(imPeak, 1) || round(gConfs(peakNum, 7)) > size(imPeak, 2)...
					|| round(gConfs(peakNum, 8)) < 1 || round(gConfs(peakNum, 7)) < 1
				% delete the configuration
				gConfs(peakNum, :) = [];
				% decrease number of found particles
				peakNum = peakNum-1;

				% jump to the end of the loop
				continue;
			end
			
			% so far warn if the center of the particle was found outside of the watershed
			% region and delete the particle
			% NOTE: this is necessary for an assumption made in activationDetection2.m where we can
			% be sure that the center of a particle can be used as an index point for the watershed region

			if imWS(locMaxY(i), locMaxX(i))~=imWS(round(gConfs(peakNum, 8)), round(gConfs(peakNum, 7)))
% 				fprintf('WARNING: fitted particle center resides outside of the corresponding watershed region, DELETING particle!\n');
% 				imWS(locMaxY(i), locMaxX(i))
% 				imWS(round(gConfs(peakNum, 8)), round(gConfs(peakNum, 7)))
% 				peakNum

				% delete the configuration
				gConfs(peakNum, :) = [];
				% decrease number of found particles
				peakNum = peakNum-1;

				% jump to the end of the loop
				continue;
			end

		end
		
		%%% do some reasonability checks

		%%% discard the peak if its amount of photons is too small
		if reasonabilityChecks(1)

			%!!!!!!!!!!!!!!!!!!!!
			% NOTE: this computation only works for my CSU setup because the photon counting is
			% specific
			%!!!!!!!!!!!!!!!!!!!!
			null = imPeak;
			% remove nan values if watershed is used
			if useWatershedTransformRegion && ~doEMFitting
				null(isnan(null)) = 0;
			end
			
			if gType == 2
				if computeNumOfParticlePhotons2D(null, dv, gConfs(peakNum, 1:1+dim-1), gConfs(peakNum, 3:3+dim-1), gConfs(peakNum, 9), photonNumInterval, bgNoiseMean, CCDSensitivity, EMGain)...
						< minNumPhotonsPerPeak
					% delete the configuration
					gConfs(peakNum, :) = [];
					% decrease number of found particles
					peakNum = peakNum-1;

					% jump to the end of the loop
					continue;
				end
			else
				if computeNumOfParticlePhotons2D(null, dv, gConfs(peakNum, 1:1+dim-1), gConfs(peakNum, 3:3+dim-1), [], photonNumInterval, bgNoiseMean, CCDSensitivity, EMGain)...
						< minNumPhotonsPerPeak
					% delete the configuration
					gConfs(peakNum, :) = [];
					% decrease number of found particles
					peakNum = peakNum-1;

					% jump to the end of the loop
					continue;
				end
			end
		end
		
		%%% check if the FWHM of the gaussian exceeded the given boundaries
		if reasonabilityChecks(2)
			if gConfs(peakNum, 3)<minMaxFWHM(1, 1)/dv(1) || gConfs(peakNum, 3)>minMaxFWHM(2, 1)/dv(1)...
					|| gConfs(peakNum, 4)<minMaxFWHM(1, 2)/dv(2) || gConfs(peakNum, 4)>minMaxFWHM(2, 2)/dv(2)
	% 			gConfs(peakNum, 3:4).*dv(1)
				% delete the configuration
				gConfs(peakNum, :) = [];
				% decrease number of found particles
				peakNum = peakNum-1;

				% jump to the end of the loop
				continue;
			end
		end
		
		%%% obtain the quality of the fitted gaussian by comparison of the fit residuals with the
		%%% noise in the image
		if reasonabilityChecks(3)
			
			fprintf('\n\nWARNING: reasonability check 3 is outdated and not used - ignoring!\n\n');
			
% 			% estimate the noise of the peak
% 			% NOTE: the noise increases with sqrt(number of photons)
% 			noise = sqrt(max(max(imPeak))./bgNoiseMean) .* (bgNoiseRMS);
% % 			noise = sqrt((gConfs(peakNum, 5)+gConfs(peakNum, 6))./bgNoiseMean) .* (bgNoiseRMS);
% % 			noise = sqrt((gConfs(peakNum, 5)+gConfs(peakNum, 6))-bgNoiseMean) .* (bgNoiseRMS);
% 			
% 			if ~useWatershedTransformRegion && ~doEMFitting
% 				% if rectangular region is used
% 				
% 				% compute the fitted gaussian
% 				fittedGaussian = gaussian2D(gConfs(peakNum,:), dv, imPeakSize(1), imPeakSize(2));
% 
% 				% compute the value of the gaussian at the position of x times the standard deviation
% 				% NOTE: this computation is actually independent of sigma
% 				% sigmaBoundaryValue = exp( -( (sigmaFactor*min(sigma)).^2 )./ (2*min(sigma)^2));
% 				sigmaBoundaryValue = exp( -( gaussianPart2Fit^2 ) / 2);
% 
% 				% duplicate the gaussian
% 				deviationMap = fittedGaussian;
% 
% 				% create binary image that has ones at the positions that are within the sigma-interval
% 				deviationMap(deviationMap<sigmaBoundaryValue) = 0;
% 				deviationMap(deviationMap>0) = 1;
% 
% 				% compute the deviation between the gaussian fit and the real data peak
% 				deviation = (imPeak-fittedGaussian);
% 
% 				% compute the mean squared error between the fits but take only an upper part of the gaussian for
% 				% comparison
% 				fitDeviation = std(deviation(deviationMap>0));
% 			else
% 				% if watershed region is used
% 
% 				% compute the fitted gaussian
% 				fittedGaussian = gaussian2D(gConfs(peakNum,:), dv, size(im, 2), size(im, 1));
% 
% 				% compute the value of the gaussian at the position of x times the standard deviation
% 				% NOTE: this computation is actually independent of sigma
% 				% sigmaBoundaryValue = exp( -( (sigmaFactor*min(sigma)).^2 )./ (2*min(sigma)^2));
% 				sigmaBoundaryValue = exp( -( gaussianPart2Fit^2 ) / 2);
% 
% 				% duplicate the gaussian
% 				deviationMap = fittedGaussian;
% 
% 				% create binary image that has ones at the positions that are within the sigma-interval
% 				deviationMap(deviationMap<sigmaBoundaryValue) = 0;
% 				deviationMap(deviationMap>0) = 1;
% 
% 				% compute the deviation between the gaussian fit and the real data peak
% 				deviation = (im-fittedGaussian);
% 
% 				% compute the mean squared error between the fits but take only an upper part of the gaussian for
% 				% comparison
% 				fitDeviation = std(deviation(deviationMap>0));
% 			end
% 			
% 			% the final fit quality is than computed from the mean squared error normalized by the noise
% 			% level
% 			fitQuality = fitDeviation/(noise*fitQualityThresholdFactor);
% 
% 			% discard a peak if it's fit is too bad, which is the case if the fit is bigger than 1
% 			% because than the mean error of the fit is higher than the noise RMS
% 			if fitQuality >1
% 				% delete the configuration
% 				gConfs(peakNum, :) = [];
% 				% decrease number of found particles
% 				peakNum = peakNum-1;
% 
% 				% jump to the end of the loop
% 				continue;
% 			end
		end

		%%% obtain the quality of the fitted gaussian by analysis of the residuals
		if reasonabilityChecks(4)
			
% 			alpha ist der Schwellwert, bei dem von h = 0 auf h = 1 gesprungen wird.
% 			Die Nullhypothese ist Normalverteilung. h = 0 bedeutet, dass die Nullhypothese einer Normalverteilung nicht abgelehnt werden kann,
% 			weil die Wahrscheinlichkeit, dass unter Annahme einer Normalverteilung solche oder ähnliche Daten entstehen zu groß (= größer als das Signifikanzlevel) ist.
% 			h = 1 bedeutet, dass die Nullhypothese einer Normalverteilung abgelehnt werden kann,
% 			weil die Wahrscheinlichkeit, dass solche oder ungewöhnlichere Daten unter Annahme einer Normalverteilung entstehen "sehr klein" ( = kleiner als Signifikanzlevel) ist.
			
			% discard a peak if it's fit is too bad, which is the case if the residuals do not
			% follow a gaussian (or maybe a poisson) distribution because than they are not just
			% varying arround the curve, no the curve is not a good fit
			% NOTE: for understanding lillitest as explained above: if we use a high threshold alpha than the probability that this data comes from a normal
			% distribution can be higher but nevertheless it is rejected as the threshold is higher
			if lillietest(residuals, fitAlphaThresh) == 1
% 			if chi2gof(deviation(deviationMap==1)) == 1
				% delete the configuration
				gConfs(peakNum, :) = [];
				% decrease number of found particles
				peakNum = peakNum-1;

				% jump to the end of the loop
				continue;
			end

% 			fprintf('\n\nWARNING: reasonability check 4 is outdated and not used - ignoring!\n\n');
% 
% 			% compute the fitted gaussian
% 			fittedGaussian = gaussian2D(gConfs(peakNum,:), dv, imPeakSize(1), imPeakSize(2));
% 
% 			% compute a map of the intensities of the gaussian that shall be used for the quality detection
% 	
% 			% compute the value of the gaussian at the position of x times the standard deviation
% 			% NOTE: this computation is actually independent of sigma
% 			% sigmaBoundaryValue = exp( -( (sigmaFactor*min(sigma)).^2 )./ (2*min(sigma)^2));
% 			sigmaBoundaryValue = exp( -( gaussianPart2Fit^2 ) / 2);
% 
% 			% duplicate the gaussian
% 			deviationMap = fittedGaussian;
% 			
% 			% create binary image that has ones at the positions that are within the sigma-interval
% 			deviationMap(deviationMap<sigmaBoundaryValue) = 0;
% 			deviationMap(deviationMap>0) = 1;
% 
% 			% compute the deviation between the gaussian fit and the real data peak
% 			deviation = (imPeak-fittedGaussian);
% 
% 			% discard a peak if it's fit is too bad, which is the case if the residuals do not
% 			% follow a gaussian (or maybe a poisson) distribution because than they are not just
% 			% varying arround the curve, no the curve is not a good fit
% % 			if lillietest(deviation(deviationMap==1)) == 1
% 			if chi2gof(deviation(deviationMap==1)) == 1
% 				% delete the configuration
% 				gConfs(peakNum, :) = [];
% 				% decrease number of found particles
% 				peakNum = peakNum-1;
% 
% 				% jump to the end of the loop
% 				continue;
% 			end
		end
		
		%%% if watershed region is used than check if the size of the minor axis of the ellipse,
		%%% covering the watershed region, has at least a size of a factor from the fwhm such that
		%%% the fit is reasonable
		if reasonabilityChecks(5)
			
			if useWatershedTransformRegion && ~doEMFitting
				
				% get size of minor axis
				null = regionprops(imWS==imWS(locMaxY(i), locMaxX(i)), 'MinorAxisLength');
				null = null.MinorAxisLength;

				% check if the length of the smallest axis is long enough
				if null < watershedRegionControlAxisSizeFactor*mean([minMaxFWHM(1, 1)/dv(1) minMaxFWHM(1, 2)/dv(2)])
					% delete the configuration
					gConfs(peakNum, :) = [];
					% decrease number of found particles
					peakNum = peakNum-1;

					% jump to the end of the loop
					continue;
				end				
			end
		end
		
% 		%%% obtain the quality of the fitted gaussian by analysis of the residuals and check against
% 		%%% the hypothesis that they follow a normal distribution
% 		if reasonabilityChecks(6)
% 			
% 			% discard a peak if it's fit is too bad, which is the case if the residuals do not
% 			% follow a normal distribution because than they are not just
% 			% varying arround the curve, no the curve is not a good fit
% % 			if lillietest(deviation(deviationMap==1)) == 1
% 			if chi2gof() == 1
% 				% delete the configuration
% 				gConfs(peakNum, :) = [];
% 				% decrease number of found particles
% 				peakNum = peakNum-1;
% 
% 				% jump to the end of the loop
% 				continue;
% 			end
% 		end
		
	end
end

% % loop over all local maximas
% for i=1:length(locMaxX)
% 	
% 	% only accept a local maxima as a peak if it is above the peak threshold
% 	if im(locMaxY(i), locMaxX(i)) > peakThresh
% 
% 		% furthermore reject peaks that are too close to the border of the image
% 		if (locMaxX(i)-imPeakSizeHalf(1))>=1 && (locMaxX(i)+imPeakSizeHalf(1))<=size(im, 2)...
% 			&& (locMaxY(i)-imPeakSizeHalf(2))>=1 && (locMaxY(i)+imPeakSizeHalf(2))<=size(im, 1)
% 
% 			% increase number of peaks
% 			peakNum = peakNum+1;
% 
% 			% get the corresponding part in the image, which is our peak where a gauss shall be fitted to
% 			imPeak = im(max(1, locMaxY(i)-imPeakSizeHalf(2)):min(size(im, 1), locMaxY(i)+imPeakSizeHalf(2))...
% 				, max(1, locMaxX(i)-imPeakSizeHalf(1)):min(size(im, 2), locMaxX(i)+imPeakSizeHalf(1)) );
% 
% 			% fit a gaussian to the peak in the image
% 			gConfs = [gConfs; fitGaussian2D(imPeak, dv, gConfInit, gConfLowerBound, gConfUpperBound, gType)];
%
% 			% update the real xy positions of the gauss center in the image
% 			gConfs(peakNum, 7) = gConfs(peakNum, 1) -1 + max(1, locMaxX(i)-imPeakSizeHalf(1));
% 			gConfs(peakNum, 8) = gConfs(peakNum, 2) -1 + max(1, locMaxY(i)-imPeakSizeHalf(2));
% 
% 			% discard the peak if its amount of photons is too small
% 			if computeNumOfParticlePhotons2D(imPeak, dv, gConfs(peakNum, 1:1+dim-1), gConfs(peakNum, 3:3+dim-1), bgNoiseMean, CCDSensitivity, EMGain)...
% 					< minNumPhotonsPerPeak
% 				% delete the configuration
% 				gConfs(peakNum, :) = [];
% 				% decrease number of found particles
% 				peakNum = peakNum-1;
%  			end
% 			
% 		end
% 	end
% 
% end

%% delete found particles that are too close together
% this might be either the same particle found twice or particles that are just too close together.
% here close is defined as less as the size of the fwhm of the brighter particle. the brighter
% particle is finally taken
% NOTE/TODO: it might happen that particles are deleted concatenated because all are close to each
% other like in a chain and thus the biggest deletes all

removeParticles = [];

% compare each particle position with all other ones and get those particles that shall be removed
for i=1:size(gConfs, 1)-1
	for j=i+1:size(gConfs, 1)
		% check which particle is the brightest and take its fwhm as the distance
		if sum(gConfs(i, 5:6), 2)>=sum(gConfs(j, 5:6), 2)
			
			if ~useConfForSPT
				% take its fwhm as the distance
				dist = gConfs(i, 3:4).*synapseDistInterval;
			else
				% take a fixed distance
				dist = [synapseDistInterval synapseDistInterval];
			end
				
			% store the index of the darker particle
			pIdx = j;
		else
			if ~useConfForSPT
				% take its fwhm as the distance
				dist = gConfs(j, 3:4).*synapseDistInterval;
			else
				% take a fixed distance
				dist = [synapseDistInterval synapseDistInterval];
			end
			
			% store the index of the darker particle
			pIdx = i;
		end
		
		% now compare the distance in all dimensions separately
		% if the distance is too low
		% NOTE: because of the discretization i do not compute the vector distance as than still the overlap of pixels would be possible
		if abs(gConfs(i, 7)-gConfs(j, 7))<dist(1) && abs(gConfs(i, 8)-gConfs(j, 8))<dist(2)
			% if both particles shall be removed
			if removeAllParticlesTooCloseToOthers == 1
				% add both indexes to the remove list
				removeParticles = [removeParticles, i, j];
			else
				% only remember the index of the darker particle
				removeParticles = [removeParticles, pIdx];
			end
		end
	end
end

%%% select those particles that remain in the list
% create full index
keepParticles = [1:size(gConfs, 1)];
% set indexes to zero that are to be removed
% if all too close particles shall be removed
if removeAllParticlesTooCloseToOthers == 1
	% remove multiple similar particle numbers before removing
	keepParticles(unique(removeParticles)) = 0;
else
	keepParticles(removeParticles) = 0;
end
% take only those indexes that are not zero and thus one has the indexes of those particles that
% remain
keepParticles = keepParticles(keepParticles~=0);

% now really remove the found particles from the list
gConfs = gConfs(keepParticles, :);

% update number of found peaks
peakNum = size(gConfs, 1);

% unload EM library
if doEMFitting
% 	unloadlibrary('pem');
% 	unloadlibrary('pem2');
end

%% draw fitting results in image
% figure, imagesc(im), colormap(gray), hold on;
% for numPeaksInitIdx = 1:size(gConfs, 1)
% 	if ~isnan(gConfs(numPeaksInitIdx, 7))
% 		% if symmetrical Gaussians shall be fitted
% 		if gType == 0
% 			drawGaussian([gConfs(numPeaksInitIdx, 7) gConfs(numPeaksInitIdx, 8) gConfs(numPeaksInitIdx, 3)/2.3548200450309 gConfs(numPeaksInitIdx, 4)/2.3548200450309 0], 'r');
% 		% else if also the angle rho can be adjusted
% 		else
% 			drawGaussian([gConfs(numPeaksInitIdx, 7) gConfs(numPeaksInitIdx, 8) gConfs(numPeaksInitIdx, 3)/2.3548200450309 gConfs(numPeaksInitIdx, 4)/2.3548200450309 gConfs(numPeaksInitIdx, 9)], 'r');
% 		end
% 	end
% end
% hold off;

%% display the results within the image
% TODO: some form of 4D representation would be necessary -> not done so far

if display
	
% 	% print results
% 	display(gConfs);

	%% mark the image at the found position of the peaks
	
	% get the original image
	imMarked = imOrig;
	
	% insert all peaks
	for i=1:peakNum

% 		%%% add a gauss at the current peak
% 		
% 		% overwrite image with zero values
% 		imMarked = zeros(size(imOrig));
% 		
% 		% compute the corresponding fitted gaussian
% 		fittedGauss = gaussian2DEllipsoid(gConfs(i,:), dv, imPeakSize(1), imPeakSize(2));
% 
% 		% add peak
% 		imMarked(round(gConfs(i, 8))-imPeakSizeHalf(2):round(gConfs(i, 8))+imPeakSizeHalf(2)...
% 			, round(gConfs(i, 7))-imPeakSizeHalf(1):round(gConfs(i, 7))+imPeakSizeHalf(1)) = fittedGauss;

		%%% add cross at the current peak
 		imMarked(round(gConfs(i, 8)), round(gConfs(i, 7))) = 0;
		imMarked(round(gConfs(i, 8))+1, round(gConfs(i, 7))) = 0;
		imMarked(round(gConfs(i, 8))-1, round(gConfs(i, 7))) = 0;
		imMarked(round(gConfs(i, 8)), round(gConfs(i, 7))+1) = 0;
		imMarked(round(gConfs(i, 8)), round(gConfs(i, 7))-1) = 0;
		
	end
	
	% insert all gaussians as gaussian shaped markers into a binary map
	imMFull = zeros(size(imOrig));
	for i=1:peakNum
		imMFull = getParticleIntensities2D(im, dv, gConfs(i, 7:7+dim-1), gConfs(i, 3:3+dim-1), 2) | imMFull;
	end
	imMFull(imMFull>0) = 1;
	imMFull(imMFull<=0) = 0;

	figure;
	overlay(im, imMFull, [0 1 0], 0.1, 'Gaussian Localization');
	

% 	% create binary map of found peaks
% 	bwSeg = zeros(size(im));
% 
% 	% load a circle image as marker
% 	if markerRadius == 2
% 		markerCircle = imread('~/data/projects/matlabCommon/circle5.png');
% 	elseif markerRadius == 3
% 		markerCircle = imread('~/data/projects/matlabCommon/circle7.png');
% 	elseif markerRadius == 4
% 		markerCircle = imread('~/data/projects/matlabCommon/circle9.png');
% 	elseif markerRadius == 5
% 		markerCircle = imread('~/data/projects/matlabCommon/circle11.png');
% 	end
% 
% 	% convert to binary
% 	markerCircle(markerCircle>0) = 1;
% 
% 	% insert all peaks
% 	for i=1:peakNum
% 
% 		% add a gauss at the current peak
% 		% compute the corresponding fitted gaussian
% 		fittedGauss = gaussian2DEllipsoid(gConfs(i,:), dv, imPeakSize(1), imPeakSize(2));
% 		
% 		% % create a segmentation mask from FWHM
% 		% fittedGauss(fittedGauss<floor(max(max(fittedGauss))/2)) = 0;
% 		% fittedGauss(fittedGauss>0) = 1;
% 		
% 		% add peak
% 		bwSeg(round(gConfs(i, 8))-imPeakSizeHalf(2):round(gConfs(i, 8))+imPeakSizeHalf(2), round(gConfs(i, 7))-imPeakSizeHalf(1):round(gConfs(i, 7))+imPeakSizeHalf(1)) = fittedGauss;
% 
% % 		% add circle at current peak
% % 		bwSeg(round(gConfs(i, 8))-markerRadius:round(gConfs(i, 8))+markerRadius, round(gConfs(i, 7))-markerRadius:round(gConfs(i, 7))+markerRadius) = markerCircle;
% 
% % 		% add point at peak
% % 		bwSeg(round(gConfs(i, 2)), round(gConfs(i, 1))) = 1;
% 
% 	end
% 
% % 	figure;
% % 	overlay(imOrig, bwSeg, [1 0 0], 0.2, 'detect gaussians');
% 	figure;
% 	surf(im);
% 	colormap(jet);
% 
% 	figure;
% 	surf(bwSeg);
% 	colormap(jet);
else
% 	% get the original image
% 	imMarked = imOrig;
	imMarked = [];
end

end