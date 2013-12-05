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


%%% try to make two very close synapses distinguishable
% usually it happens that two close synapses are fitted as one together. this shall be avoided via
% creating areas arround the local maxima (using watershed) such that each synapse is taken as its
% own
if useWatershedTransformRegion
	% invert the image first so that minima become maxima and vice versa
% 	imWS = cast(abs(im-(max(max(im)))), 'uint16');
	imWS = cast(max(max(im))-im, 'uint16');
	% do watersheed
	imWS = watershed(imWS, 8);

% 	% plot watershed result
% 	figure, imshow(label2rgb(imWS));
end

%%% detect possible peaks that will be fitted later on

% find local maximas in the image and save as binary
if isempty(useImTopHatFilter)
	% find all local maxima
	imLocMax = imregionalmax(im, 8);
else
% 	stdSignalValues = std(signalValues)
% 	peakThresh
	% find local maxima using top hat filtering
	imLocMax = topHatFilter(im, std(signalValues).*useImTopHatFilter(1), useImTopHatFilter(2), useImTopHatFilter(3), 1 , 1);
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
