% unsaved changes:
% -now all frames where no particles could be found are really put together in one track
% -the diagonal length of a voxel is now taken for maxFrameParticleDistance
% -now even leftover tracks are tried to be connected between each other
%
% function trackParticles2D(im, preprocess, display, EMGain, frameRate, speed)
%
% tracks particles over several frames
% NOTE: thsi tracking algorithms assumes that tracks don't overlap thus it optimizes for long tracks
% by first tacking complete ones, than filling up all blinking ones and finally even leftover tracks
% are tried to fit to together to longer tracks
%
% parameters:
% -im			- the images containing the particles
% -dv			- real size in each dimension of one voxel
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
% 		[08] - gaussianPart2Fit -> 1/3
%		[09] - intensityOutliers -> 0.1
%		[10] - reasonabilityChecks -> [1 1 1]
%		[11] - subImSizefactor -> 2
%		[12] - useWatershedTransformRegions -> 0
%		[13] - weightedLsqPoisson -> [] (empty matrix means don't use weighted LSQ, put a number for adjusting, so far the number doesn't matter)
%		trackingParams{1}=3,trackingParams{2}=1,trackingParams{3}=[200, 200; 1400 1400],trackingParams{4}=1,trackingPara
%		ms{5}=2,trackingParams{6}=25,trackingParams{7}=2,trackingParams{8}=1/3,trackingParams{9}=0.1
%		, trackingParams{10}=[1 1 1], trackingParams{11}=2, trackingParams{12}=1,
%		trackingParams{13}=[]
% -preprocess	- 1 if the image shall be preprocessed, 0 otherwise
% -display		- 1 if debug and result information shall be created, 0 otherwise
% -EMGain		- EMGain used during acqisition
% -frameRate	- frame rate of the image aquisition
% %-speed		- maximal speed of the particle to track (the speed hat to be provided in nm/s (usually 40000))
% -binning		- level of binning used
%
% returns:
% -imMarked		-if display is turned on the last frame will
function [imMarked results] = trackParticles2DBetha(im, dv, trackingParams, preprocess, display, EMGain, frameRate, binning)

%% constructor

% compute number of dimensions (just that i don't have to use the number all time)
dim = length(dv);
% correct the voxel size if binning occurs
dv(1)=dv(1)*binning;
dv(2)=dv(2)*binning;
% compute number of frames
maxFrameNum = size(im, dim+1);
if display
	% allocate memory for the marked output image
	imMarked = cast(zeros(size(im)), 'uint16');
end

%% parameters/variables

% if results shall be displayed
printResults = display;

% -plots		- defines what debugging options shall be provided:
% 			(1) - print fitted qds over image
plots = [0];

% size of the region to plot the fitted qds over the image as a factor of the standard deviation
plottingStdFactor = 2;

% number of photo-electrons per count of the A/D converter
CCDSensitivity = 12.4;

% number of voxels that are valid to search for a successor particle in the next frame
maxVoxelDisplacement = trackingParams{1};

% maximal distance consecutive particles may have between frames
% NOTE: the value is taken as the euclidean distance WITHOUT extracting the sqrt() to save
% performance
% NOTE: the value has to be in real size
% TODO: make this less arbitrary
% maxFrameParticleDist = sum(dv.^2);	% so far its simply the euclidean distance of one voxel
% maxFrameParticleDist = (speed/(frameRate/size(im, dim)))^2; % 200nm^2 per second is the usual value for neurexin
% maxFrameParticleDist = (maxVoxelDisplacement*dv(1)/binning)^2	% simply x times the length of the x direction
% maxFrameParticleDist = sum((maxVoxelDisplacement.*dv./binning).^2)./dim;	% mean real size of a voxel
maxFrameParticleDist = sum((maxVoxelDisplacement.*dv./binning).^2);	% mean real diagonal of a voxel

% factor for the size of the synapse within which no other synapse is allowed to occure
synapseDistInterval = trackingParams{2};

% minimal and maximum values in nm for the FWHM to be accepted for a gaussian fit. example: [minx miny;maxx maxy]
minMaxFWHM = trackingParams{3};

% type of the gaussian:
% -gType		- (0) -> gaussian with possibly varying width in x and y direction
%				- (1) -> gaussian with the same width in x and y direction
%				- (2) -> elliptical gaussian having varying x and y widths but can also be turned by
%				an angle theta
gType = trackingParams{4};

% list containing all found particles over all frames
particles = [];
% list containing the accuracies for all found particles over all frames
accuracy = [];
% list containing the SNRs for all frames (this includes the foreground and background version)
snr = [];
% list containing the mean accuracy for each frame
meanAccuracy = [];
% index connecting frames to all particles within that frame
frame2particleIdx = {};
% index connecting all particles to their frame
particle2frameIdx = {};
% index connecting all tracks to their corresponding particles (these are the partial tracks like they come
% out of the pure tracking procedure)
partialTrack2particleIdx = {};
% complete track to particle index (these are those tracks that reach from the first to the last
% frame without breaks)
completeTrack2particleIdx = {};
% blinking track to particle index (these tracks start form the beginning but blink and don't need
% to last to the last frame)
blinkingTrack2particleIdx = {};
% leftover track to particle index (these are parts of tracks that are unused and unconnected so far)
leftoverTrack2particleIdx = {};

% results as cell array
results = {};

%% preliminaries

%% main

%% track particles over frames

% iterate over frames and create list of particles over all frames and create the according indexes
for i=1:maxFrameNum
	
	fprintf('[%d] ', i);
	
	% detect all particles for a frame
	[fParticles imMarked(:,:,i)] = locateParticles2DBetha(im(:,:,i), dv, trackingParams, preprocess, display, CCDSensitivity, EMGain);

	% if no particles were found
	if isempty(fParticles)
		fParticles = [-1, -1, -1, -1, -1, -1, -1, -1];
		faccuracy = [-1 -1 -1 -1];
		fmeanAcc = [-1 -1 -1 -1];
		fminAcc = [-1 -1 -1 -1];
		fmaxAcc = [-1 -1 -1 -1];
		ffgSNR = -1;
		fbgSNR = -1;
	else
		%!!!!!!!!!!!!!!!!!!!!
		% NOTE: these computations only work for my CSU setup because the photon counting is
		% specific
		%!!!!!!!!!!!!!!!!!!!!

		% compute statistics for the accuracy (in real size) of the particles within the stack
		if gType == 2
			[faccuracy fmeanAcc fminAcc fmaxAcc] = computeLocalizationAccuracy2D(im(:,:,i), dv...
				, fParticles(:, 7:7+dim-1), fParticles(:, 3:3+dim-1), fParticles(:, 9), CCDSensitivity, EMGain);
		else
			[faccuracy fmeanAcc fminAcc fmaxAcc] = computeLocalizationAccuracy2D(im(:,:,i), dv...
				, fParticles(:, 7:7+dim-1), fParticles(:, 3:3+dim-1), [], CCDSensitivity, EMGain);
		end
		
		if isnan(fmeanAcc)
			fmeanAcc = [-1 -1 -1 -1];
			fminAcc = [-1 -1 -1 -1];
			fmaxAcc = [-1 -1 -1 -1];
		end
		
		% estimate signal to noise ratio
		% 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		% ich glaube bei cheezum2002, kubitscheck2000 und rogers2006 nehmen die einfach die mittlere
		% (max)peak pixel intensitaet als intensitaet und dessen standard deviation als rauschen und cih
		% bisher ja die mittlere intensitaet wie cheezum2002 das benennt (muss ich mich einfach zu einem entscheiden)
		% ZUDEM sollte ich ueberlegen ob der algorithmus so richtig ist weil ich ja quasie das rauschen
		% ueber die partikel hinweg mit messe und letztendlich sind QDs nunmal nicht immer gleich hell
		% sprich mein SNR ist zu niedrig berechnet da das rauschen zu hoch geschaetzt wird. AM BESTEN wenns
		% gebraucht wird NEU SCHREIBEN!!!!!!!!
		% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if gType == 2
			[ffgSNR fbgSNR] = estimateSNR2D(cast(round(im(:,:,i)), 'uint16'), dv, fParticles(:, 7:7+dim-1), fParticles(:, 3:3+dim-1), fParticles(:, 5), fParticles(:, 6), fParticles(:, 9));
		else
			[ffgSNR fbgSNR] = estimateSNR2D(cast(round(im(:,:,i)), 'uint16'), dv, fParticles(:, 7:7+dim-1), fParticles(:, 3:3+dim-1), fParticles(:, 5), fParticles(:, 6), []);
		end
		
% 		if isnan(ffgSNR)
% 			ffgSNR = -1;
% 			fbgSNR = -1;
% 		end
		for k = 1:size(ffgSNR, 1)
			if isnan(ffgSNR(k))
				ffgSNR(k) = -1;
				fbgSNR(k) = -1;
			end
		end
	end

	% create full list of particles
	particles = [particles; fParticles];

	% create full list of accuracies
	accuracy = [accuracy; faccuracy];

	% create list of mean accuracies per frame
	meanAccuracy = [meanAccuracy; fmeanAcc];
	
	% create full list of SNRs
	snr =  [snr; [ffgSNR fbgSNR]];
	
	% create frame to particle index
	if i == 1
		frame2particleIdx(i) = {[1:size(fParticles, 1)]};
	else
		frame2particleIdx(i) = {[1:size(fParticles, 1)]+frame2particleIdx{i-1}(end)};
	end

	% create particle to frame index
	% NOTE: i use a cell array as well even that it is unnecessary because of simmilar usage of all
	% indexes
	for j = 1:size(fParticles, 1)
		particle2frameIdx = [particle2frameIdx, i];
	end
	
	% if we want to plot the fitting results over the image
	if plots(1) && i < 4
		
		% compute mask covering the fitted qds
		IMask = zeros(size(im, 1), size(im, 2));
		
		% loop over particles and set the according number at pixels of the same region in the mask
		for pNum = 1 : size(fParticles, 1)

			% compute gaussian
			if gType == 2
				null = gaussian2DEllipse([fParticles(pNum, 7), fParticles(pNum, 8), fParticles(pNum, 3), fParticles(pNum, 4), 1, 0, 0, 0, fParticles(pNum, 9)], dv, size(im, 2), size(im, 1));
			else
				null = gaussian2D([fParticles(pNum, 7), fParticles(pNum, 8), fParticles(pNum, 3), fParticles(pNum, 4), 1, 0], dv, size(im, 2), size(im, 1));
			end

			% cut the gaussian at a certain sigma-interval and insert the pNum as the index for the
			% region
			% NOTE: this computation is actually independent of sigma
			null(null<exp( -( plottingStdFactor^2 ) / 2)) = 0;
			null(null>0) = pNum;

			% add to the intensity mask
			IMask = IMask + null;

			% test that each region was added only once
			if sum(sum(IMask==pNum))~=sum(sum(null==pNum))
				fprintf('WARNING: two particles felt into the same region, ignoring this so far\n');
			end
		end

		% also create a binary intensity mask that has simply ones at all particle positions
		IMaskBinary = IMask;
		IMaskBinary(IMaskBinary>0) = 1;

		% mask on qd image
		figure;
		overlay(imresize(im(:,:,i), 1.5, 'nearest'), imresize(IMaskBinary, 1.5, 'nearest'), [0 1 0], 0.2, sprintf('gaussian mask on qd image depth: %d', i));
	end

end
fprintf('\n');

%% create partial tracks
%% perform tracking thus connect particles over frames to tracks. the result is a track to
%% particles index
% NOTE: the tracking index/algorithm is simple so far:
% -it is assumed that particles are immobilised and thus simply the closest particle in the next frame is taken.
% this is the first particle that matches.
% -if a track stops because of blinking or just because the particle has not been found in the next frame
% than this track is closed for ever.
% -new tracks at later frames have the index of the head of the track for those frames before.
% -closed tracks always continue to have the index of the tail for all remaining frames
% one can distinguish just started tracks and closed tracks by checking if the last particle index
% in teh track really belongs to that frame
for i=1:maxFrameNum

	% for the first frame
	if i == 1
		% just take the particles of the first frame as the start of the tracks
		for j=1:length(frame2particleIdx{1})
			partialTrack2particleIdx(j) = {frame2particleIdx{1}(j)};
		end
		
	% for all consecutive frames
	else
		% get list of all particles that are the tail of the current tracks
		pre = [];

		for j=1:size(partialTrack2particleIdx, 2)
			pre = [pre; particles(partialTrack2particleIdx{j}(end), :)];
		end

		% get list of all particles of the current frame
		curr = particles(frame2particleIdx{i}, :);
		
		%%% loop over preceeding particles and try to find a corresponding particle in the following
		%%% frame by first finding the closest one and than deciding if it is close enough
		for j=1:size(pre, 1)
			
			% if the track is not already closed
			% NOTE: this is indicated by simmilar particle indexes in consecutive frames thus
			% if its just the second frame as well as if the two last track indexes are unequal as well as if the last two track indexes are
			% unequal but the last index is a correct index for the frame (this is required to
			% distinguish between closed and just started tracks)
			if (i == 2) ||...
				(partialTrack2particleIdx{j}(end)~=partialTrack2particleIdx{j}(end-1)) ||...
				((partialTrack2particleIdx{j}(end)==partialTrack2particleIdx{j}(end-1)) && (particle2frameIdx{partialTrack2particleIdx{j}(end)}==i-1))
			
				% reset distance vector
				minDist = [];

				% loop over current particles
				for k=1:size(curr, 1)

					% compute current distance in all dimensions
					% NOTE: taking the voxel size into account is important here if they are not
					% uniquely spaced
					dimDist = abs((pre(j, 7:7+dim-1).*dv)-(curr(k, 7:7+dim-1).*dv));

					% compute current euclidean distance and store it in a distance vector
					% NOTE: to avoid unnecessary sqrt computation we compare the distances without
					% sqrt()
					minDist(k) = sum(dimDist.^2);
				end

				% take the distance and the index of the closest particle
				[minDist closestParticle] = min(minDist);

				% check if the particle is close enough to be a valid successor
				if minDist<=maxFrameParticleDist
					
					% than it is added to the track to particle index
					partialTrack2particleIdx(j) = {[partialTrack2particleIdx{j}, frame2particleIdx{i}(closestParticle)]};

					% and deleted from the list of available particals in the current frame
					% NOTE: "deletion" means that it gets a huge negative distance assigned such that it
					% falls out as a possible candidate (this is at least -2)
					curr(closestParticle, 7:7+dim-1) = -2*size(im(:,:,i));

				else
					% this particular track is closed and marked with the same particle index as the
					% last one
					partialTrack2particleIdx(j) = {[partialTrack2particleIdx{j}, partialTrack2particleIdx{j}(end)]};
				end

			% if this track was already closed	
			else
				% just add the particle index of the last position of the track again
				partialTrack2particleIdx(j) = {[partialTrack2particleIdx{j}, partialTrack2particleIdx{j}(end)]};
			end
		end
		
		%%% by now every track is finished by either adding a consecutive particle or closing the track
		%%% now those current particles that couldn't be added to a track make up a new track
		%%% theirself

		% loop over current particles and add those that have non-negative positions as new tracks
		% NOTE: -1 ones which means that none was found are added nevertheless
		for k=1:size(curr, 1)
			if curr(k, 7)>=-1
				% create new track with all values having the same particle index as the current one
				% that starts the track
				partialTrack2particleIdx = [partialTrack2particleIdx, zeros(1, i)+frame2particleIdx{i}(k)];
			end
		end

	end
end

% for null = 1:length(partialTrack2particleIdx)
% 	partialTrack2particleIdx{null}
% end

%% create full tracks from partial ones
% NOTE: thus blinking is treated
% NOTE: this algorithm is simple so far:
% -blinking tracks are so far just bild from the distance criterium. if the particle would move it
% would moust probably fail to build a consistent blinking track
% -only tracks that start from the beginning are build up to complete or blinking tracks. partial
% tracks inbetween are not connected so far

% copy the partialTrack2particleIdx
leftoverTrack2particleIdx = partialTrack2particleIdx;

%%% first collect all full tracks

% loop over all leftover tracks (from behind)
for i = fliplr(1:size(leftoverTrack2particleIdx, 2))
	% if there are only unique particle indexes in the track than it is a full one and stored in the
	% complete track list
	if length(leftoverTrack2particleIdx{i}) == length(unique(leftoverTrack2particleIdx{i}))
		% store it in completed track list
		completeTrack2particleIdx = [completeTrack2particleIdx, leftoverTrack2particleIdx{i}];
		% delete it from the leftover list
		leftoverTrack2particleIdx(i) = [];
	end
end

%%% now create list of blinking tracks that really start from the beginning

% loop over remaining leftover tracks (from behind)
for i = fliplr(1:size(leftoverTrack2particleIdx, 2))
	% if the first particle index really belongs to the first frame than this is a blinking track that
	% starts from the beginning
	if particle2frameIdx{leftoverTrack2particleIdx{i}(1)} == 1
		% store it in blinking track list
		blinkingTrack2particleIdx = [blinkingTrack2particleIdx, leftoverTrack2particleIdx{i}];
		% delete it from the leftover list
		leftoverTrack2particleIdx(i) = [];
	end
end

%% try to find leftover tracks that might belong to one of the started blinking tracks

% % first sort the blinking tracks NOTE: that doesn't definately make sense because maybe its better
% % to take long lists first such that the blinking phases become short
% blinkingTrack2particleIdx = sortBlinkingTrackLists(blinkingTrack2particleIdx, particle2frameIdx);

% thus loop over blinking tracks
for i = 1:size(blinkingTrack2particleIdx, 2)

	% while the track is not fully finished try to find further consecutive tracks that belong to
	% this one
	% not fully finished means that the last particle index is not really a particle of the last
	% frame in the image
	% NOTE: if no consecutive tracks are found the while loop will be left using a break statement
	while particle2frameIdx{blinkingTrack2particleIdx{i}(end)} ~= maxFrameNum
	
		% get the last correct frame after which the track got lost
		lastFrame = particle2frameIdx{blinkingTrack2particleIdx{i}(end)};

		% reset distance vector
		minDist = [];

		% loop over remaining leftover list and try to find a valid successor of the started blinking
		% track via a distance criterion
		for j = 1:size(leftoverTrack2particleIdx, 2)

			% if the leftover track starts with a particle index not earlier and the same than the last frame of the
			% blinking track
			if particle2frameIdx{leftoverTrack2particleIdx{j}(1)} > lastFrame
				% compute current distance in all dimensions
				% NOTE: taking the voxel size into account is important here if they are not
				% uniquely spaced
				dimDist = abs((particles(leftoverTrack2particleIdx{j}(1), 7:7+dim-1).*dv)...
					-(particles(blinkingTrack2particleIdx{i}(lastFrame), 7:7+dim-1).*dv));
			% if the leftover track starts too early
			else
				% just assign a huge distance in all dimensions such that this particle is sorted out
				% by the distance criterion
				dimDist = 2*(size(im(:,:,1)).*dv);
			end

			% compute current euclidean distance and store it in a distance vector
			% NOTE: to avoid unnecessary sqrt computation we compare the distances without
			% sqrt()
			minDist(j) = sum(dimDist.^2);
		end

		% detect particles that are within the valid distance
		closeParticles = minDist-maxFrameParticleDist;
		closeParticles(closeParticles<=0) = -1;
		closeParticles(closeParticles>0) = 0;
				
		% if more than one particle is within the valid distance
		if abs(sum(closeParticles))>1
			% reset lowest frame index
			lowestFrame = maxFrameNum+1;

			% choose the one to be the successor that restarts at the lowest frame
			
			% thus loop over close particles
			for k = find(closeParticles)
				% and search for the lowest one
				if lowestFrame > particle2frameIdx{leftoverTrack2particleIdx{k}(1)}
					% remember current lowest frame
					lowestFrame = particle2frameIdx{leftoverTrack2particleIdx{k}(1)};
					% as well as the regarding index in the leftover list
					closestParticle = k;
				end
			end
			
			% then this track is appended to the current blinking track
			blinkingTrack2particleIdx(i) = {[blinkingTrack2particleIdx{i}(1:particle2frameIdx{leftoverTrack2particleIdx{closestParticle}(1)}-1)...
				, leftoverTrack2particleIdx{closestParticle}(particle2frameIdx{leftoverTrack2particleIdx{closestParticle}(1)}:end)]};

			% and deleted from the list of leftover tracks
			leftoverTrack2particleIdx(closestParticle) = [];
			
		% if at most one is within the valid distance
		else
			% take the distance and the index of the closest particle
			[minDist closestParticle] = min(minDist);

			% check if the particle is close enough to be a valid successor
			if minDist<=maxFrameParticleDist

				% then this track is appended to the current blinking track
				blinkingTrack2particleIdx(i) = {[blinkingTrack2particleIdx{i}(1:particle2frameIdx{leftoverTrack2particleIdx{closestParticle}(1)}-1)...
					, leftoverTrack2particleIdx{closestParticle}(particle2frameIdx{leftoverTrack2particleIdx{closestParticle}(1)}:end)]};

				% and deleted from the list of leftover tracks
				leftoverTrack2particleIdx(closestParticle) = [];
			% if no consecutive track was found than leave the while loop
			else
				break; % leave while loop
			end
		end
	end
end

%% try to connect leftover tracks between each other

% copy the leftovertrack2particle idx and delete it
restTrack2particleIdx = leftoverTrack2particleIdx;
leftoverTrack2particleIdx = {};

% set index for leftovertrack that will be build up from the rest track to start value
k = 0;

% as long as their are tracks to connect
while ~isempty(restTrack2particleIdx)

	% increase leftover track index
	k = k+1;
	
	% reset earliest resttrack index
	earliestRestTrack = 0;
	
	% take the list from the rest list that starts first
	for i = 1:length(restTrack2particleIdx)
		
		if i == 1
			% take the first anyway
			earliestRestTrack = i;
		else
			% ckeck if the current earliest rest track starts later than the next rest track
			if particle2frameIdx{restTrack2particleIdx{earliestRestTrack}(1)} > particle2frameIdx{restTrack2particleIdx{i}(1)}
				% if so than remember the index of the next rest track
				earliestRestTrack = i;
			end
		end
	end
	
	% copy the earliest resttrack into the leftover track and delete it from the rest list
	leftoverTrack2particleIdx(k) = restTrack2particleIdx(earliestRestTrack);
	restTrack2particleIdx(earliestRestTrack) = [];
		
	%%% while the current leftover track is not fully finished try to find further consecutive rest tracks
	% that belong to this one
	% not fully finished means that the last particle index is not really a particle of the last
	% frame in the image
	% NOTE: if no consecutive tracks are found the while loop will be left using a break statement
	while particle2frameIdx{leftoverTrack2particleIdx{k}(end)} ~= maxFrameNum
	
		% get the last correct frame after which the track got lost
		lastFrame = particle2frameIdx{leftoverTrack2particleIdx{k}(end)};

		% reset distance vector
		minDist = [];

		% loop over remaining rest list and try to find a valid successor of the current leftover
		% track via a distance criterion
		for j = 1:size(restTrack2particleIdx, 2)

			% if the rest track starts with a particle index not earlier and the same than the last frame of the
			% leftover track
			if particle2frameIdx{restTrack2particleIdx{j}(1)} > lastFrame
				% compute current distance in all dimensions
				% NOTE: taking the voxel size into account is important here if they are not
				% uniquely spaced
				dimDist = abs((particles(restTrack2particleIdx{j}(1), 7:7+dim-1).*dv)...
					-(particles(leftoverTrack2particleIdx{k}(lastFrame), 7:7+dim-1).*dv));
			% if the rest track starts too early
			else
				% just assign a huge distance in all dimensions such that this particle is sorted out
				% by the distance criterion
				dimDist = 2*(size(im(:,:,1)).*dv);
			end

			% compute current euclidean distance and store it in a distance vector
			% NOTE: to avoid unnecessary sqrt computation we compare the distances without
			% sqrt()
			minDist(j) = sum(dimDist.^2);
		end

		% detect particles that are within the valid distance
		closeParticles = minDist-maxFrameParticleDist;
		closeParticles(closeParticles<=0) = -1;
		closeParticles(closeParticles>0) = 0;
				
		% if more than one particle is within the valid distance
		if abs(sum(closeParticles))>1
			% reset lowest frame index
			lowestFrame = maxFrameNum+1;

			% choose the one to be the successor that restarts at the lowest frame
			
			% thus loop over close particles
			for i = find(closeParticles)
				% and search for the lowest one
				if lowestFrame > particle2frameIdx{restTrack2particleIdx{i}(1)}
					% remember current lowest frame
					lowestFrame = particle2frameIdx{restTrack2particleIdx{i}(1)};
					% as well as the regarding index in the leftover list
					closestParticle = i;
				end
			end
			
			% then this track is appended to the current leftover track
			leftoverTrack2particleIdx(k) = {[leftoverTrack2particleIdx{k}(1:particle2frameIdx{restTrack2particleIdx{closestParticle}(1)}-1)...
				, restTrack2particleIdx{closestParticle}(particle2frameIdx{restTrack2particleIdx{closestParticle}(1)}:end)]};

			% and deleted from the list of rest tracks
			restTrack2particleIdx(closestParticle) = [];
			
		% if at most one is within the valid distance
		else
			% take the distance and the index of the closest particle
			[minDist closestParticle] = min(minDist);

			% check if the particle is close enough to be a valid successor
			if minDist<=maxFrameParticleDist

				% then this track is appended to the current leftover track
				leftoverTrack2particleIdx(k) = {[leftoverTrack2particleIdx{k}(1:particle2frameIdx{restTrack2particleIdx{closestParticle}(1)}-1)...
					, restTrack2particleIdx{closestParticle}(particle2frameIdx{restTrack2particleIdx{closestParticle}(1)}:end)]};

				% and deleted from the list of rest tracks
				restTrack2particleIdx(closestParticle) = [];
			% if no consecutive track was found than leave the while loop
			else
				break; % leave while loop
			end
		end
	end
end

% for null = 1:length(leftoverTrack2particleIdx)
% 	leftoverTrack2particleIdx{null}
% end

% disp('ACCCURACY EVALUATION:');

completeTrackMean = [];
for i=1:size(completeTrack2particleIdx, 2)
	completeTrackMean = [completeTrackMean; mean(particles(completeTrack2particleIdx{i}, 7:8), 1)];
	if printResults == 1
		fprintf('complete position means [%d/%d]: %fpixel  %fpixel\n', i, length(unique(completeTrack2particleIdx{i})), completeTrackMean(i, 1), completeTrackMean(i, 2));
	end
end
blinkingTrackMean = [];
for i=1:size(blinkingTrack2particleIdx, 2)
	blinkingTrackMean = [blinkingTrackMean; mean(particles(unique(blinkingTrack2particleIdx{i}), 7:8), 1)];
	if printResults == 1
		fprintf('blinking position means [%d/%d]: %fpixel  %fpixel\n', i, length(unique(blinkingTrack2particleIdx{i})), blinkingTrackMean(i, 1), blinkingTrackMean(i, 2));
	end
end
leftoverTrackMean = [];
for i=1:size(leftoverTrack2particleIdx, 2)
	leftoverTrackMean = [leftoverTrackMean; mean(particles(unique(leftoverTrack2particleIdx{i}), 7:8), 1)];
 	if printResults == 1
		fprintf('leftover position means [%d/%d]: %fpixel  %fpixel\n', i, length(unique(leftoverTrack2particleIdx{i})), leftoverTrackMean(i, 1), leftoverTrackMean(i, 2));
	end
end

completeTrackSTD = [];
for i=1:size(completeTrack2particleIdx, 2)
	completeTrackSTD = [completeTrackSTD; std(particles(completeTrack2particleIdx{i}, 7:8), 0, 1).*dv];
 	if printResults == 1
		fprintf('complete standard deviation [%d/%d]: %fnm  %fnm\n', i, length(unique(completeTrack2particleIdx{i})), completeTrackSTD(i, 1), completeTrackSTD(i, 2));
	end
end
blinkingTrackSTD = [];
for i=1:size(blinkingTrack2particleIdx, 2)
	blinkingTrackSTD = [blinkingTrackSTD; std(particles(unique(blinkingTrack2particleIdx{i}), 7:8), 0, 1).*dv];
 	if printResults == 1
		fprintf('blinking standard deviation [%d/%d]: %fnm  %fnm\n', i, length(unique(blinkingTrack2particleIdx{i})), blinkingTrackSTD(i, 1), blinkingTrackSTD(i, 2));
	end
end
leftoverTrackSTD = [];
for i=1:size(leftoverTrack2particleIdx, 2)
	leftoverTrackSTD = [leftoverTrackSTD; std(particles(unique(leftoverTrack2particleIdx{i}), 7:8), 0, 1).*dv];
 	if printResults == 1
		fprintf('leftover standard deviation [%d/%d]: %fnm  %fnm\n', i, length(unique(leftoverTrack2particleIdx{i})), leftoverTrackSTD(i, 1), leftoverTrackSTD(i, 2));
	end
end

completeTrackFWHM = [];
for i=1:size(completeTrack2particleIdx, 2)
	completeTrackFWHM = [completeTrackFWHM; mean(particles(completeTrack2particleIdx{i}, 3:4), 1).*dv];
 	if printResults == 1
		fprintf('complete FWHMs [%d/%d]: %fnm  %fnm\n', i, length(unique(completeTrack2particleIdx{i})), completeTrackFWHM(i, 1), completeTrackFWHM(i, 2));
	end
end
blinkingTrackFWHM = [];
for i=1:size(blinkingTrack2particleIdx, 2)
	blinkingTrackFWHM = [blinkingTrackFWHM; mean(particles(unique(blinkingTrack2particleIdx{i}), 3:4), 1).*dv];
 	if printResults == 1
		fprintf('blinking FWHMs [%d/%d]: %fnm  %fnm\n', i, length(unique(blinkingTrack2particleIdx{i})), blinkingTrackFWHM(i, 1), blinkingTrackFWHM(i, 2));
	end
end
leftoverTrackFWHM = [];
for i=1:size(leftoverTrack2particleIdx, 2)
	leftoverTrackFWHM = [leftoverTrackFWHM; mean(particles(unique(leftoverTrack2particleIdx{i}), 3:4), 1).*dv];
	if printResults == 1
		fprintf('leftover FWHMs [%d/%d]: %fnm  %fnm\n', i, length(unique(leftoverTrack2particleIdx{i})), leftoverTrackFWHM(i, 1), leftoverTrackFWHM(i, 2));
	end
end

completeTrackThompson = [];
for i=1:size(completeTrack2particleIdx, 2)
	completeTrackThompson = [completeTrackThompson; mean(accuracy(completeTrack2particleIdx{i}, 1:4), 1)];
 	if printResults == 1
		fprintf('complete mean thompson accuracy [%d/%d]: (simple) %fnm  %fnm    (soph.) %fnm  %fnm\n'...
 		, i, length(unique(completeTrack2particleIdx{i})), completeTrackThompson(i, 1), completeTrackThompson(i, 2), completeTrackThompson(i, 3), completeTrackThompson(i, 4));
	end
end
blinkingTrackThompson = [];
for i=1:size(blinkingTrack2particleIdx, 2)
	blinkingTrackThompson = [blinkingTrackThompson; mean(accuracy(unique(blinkingTrack2particleIdx{i}), 1:4), 1)];
 	if printResults == 1
		fprintf('blinking mean thompson accuracy [%d/%d]: (simple) %fnm  %fnm    (soph.) %fnm  %fnmm\n'...
 		, i, length(unique(blinkingTrack2particleIdx{i})), blinkingTrackThompson(i, 1), blinkingTrackThompson(i, 2), blinkingTrackThompson(i, 3), blinkingTrackThompson(i, 4));
	end
end
leftoverTrackThompson = [];
for i=1:size(leftoverTrack2particleIdx, 2)
	leftoverTrackThompson = [leftoverTrackThompson; mean(accuracy(unique(leftoverTrack2particleIdx{i}), 1:4), 1)];
	if printResults == 1
		fprintf('leftover mean thompson accuracy [%d/%d]: (simple) %fnm  %fnm    (soph.) %fnm  %fnm\n'...
		, i, length(unique(leftoverTrack2particleIdx{i})), leftoverTrackThompson(i, 1), leftoverTrackThompson(i, 2), leftoverTrackThompson(i, 3), leftoverTrackThompson(i, 4));
	end
end

% delete empty frames
prepAccuracy = accuracy;
for i=fliplr(1:size(prepAccuracy, 1))
	if prepAccuracy(i, 1) == -1
		prepAccuracy(i,:) = [];
	end
end
meanThompson = mean(prepAccuracy, 1);
 if printResults == 1
		fprintf('mean thompson accuracy: (simple) %fnm    %fnm    (soph.) %fnm    %fnm\n'...
 	, meanThompson(1), meanThompson(2), meanThompson(3), meanThompson(4));

 	fprintf('simple: minX: %fnm  maxX: %fnm    minY: %fnm  maxY: %fnm \n'...
 		, min(prepAccuracy(:,1)), max(prepAccuracy(:,1)), min(prepAccuracy(:,2)), max(prepAccuracy(:,2)));
end

% delete empty frames
prepSnr = snr;
for i=fliplr(1:size(prepSnr, 1))
	if prepSnr(i, 1) == -1
		prepSnr(i,:) = [];
	end
end
meanSNR = mean(prepSnr, 1);
if printResults == 1
		fprintf('mean snr: (fg) %f    (bg) %f\n', meanSNR(1), meanSNR(2));

 	fprintf('(fg): minSNR: %f    maxSNR: %fnm    (bg): minSNR: %f    maxSNR: %fnm\n', min(prepSnr(:,1)), max(prepSnr(:,1)), min(prepSnr(:,2)), max(prepSnr(:,2)));
end


% disp('particles:');
% particles(:, 7)=particles(:, 7)/100;
% particles
% disp('accuracies:');
% accuracy
% disp('partial tracks:');
% partialTrack2particleIdx{:}
% disp('complete tracks:');
% completeTrack2particleIdx{:}
% disp('blinking tracks:');
% blinkingTrack2particleIdx{:}
% disp('leftover tracks:');
% leftoverTrack2particleIdx{:}

% track = particles(unique(blinkingTrack2particleIdx{1}),9:11);
% plot3(track(:, 1)*dv(1), track(:, 2)*dv(2), track(:, 3)*dv(3));

% create results vector
results(1) = {particles};
results(2) = {accuracy};
results(3) = {snr};
results(4) = {meanAccuracy};
results(5) = {frame2particleIdx};
results(6) = {particle2frameIdx};
results(7) = {partialTrack2particleIdx};
results(8) = {completeTrack2particleIdx};
results(9) = {blinkingTrack2particleIdx};
results(10) = {leftoverTrack2particleIdx};
results(11) = {completeTrackMean};
results(12) = {blinkingTrackMean};
results(13) = {leftoverTrackMean};
results(14) = {completeTrackSTD};
results(15) = {blinkingTrackSTD};
results(16) = {leftoverTrackSTD};
results(17) = {completeTrackFWHM};
results(18) = {blinkingTrackFWHM};
results(19) = {leftoverTrackFWHM};
results(20) = {completeTrackThompson};
results(21) = {blinkingTrackThompson};
results(22) = {leftoverTrackThompson};
results(23) = {meanThompson};
results(24) = {meanSNR};

end

% subfunction list = sortBlinkingTrackLists(list)
%
% sorts the lists such that the shortest one comes first using simple bubblesort
function list = sortBlinkingTrackLists(list, particle2frame)

% loop over all tracks
for i = 1:size(list, 2)-1
	for j = i:size(list, 2)-1
		if particle2frame{list{j}(end)} > particle2frame{list{j+1}(end)}
			tmp = list{j+1};
			list(j+1) = {list{j}};
			list(j) = {tmp};
		end
	end
end

end
