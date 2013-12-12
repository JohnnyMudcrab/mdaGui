% parameters:
%	im - 3d image
%	dv - real dimensions of a voxel
%	trackingParams - tracking parameters
%	plots - set it to 0
%	EMGain - set it to 64
%	p - cell array containing gConfs for the particles in each frame
%	fwhmConfX - parameters for the fwhm curve in x direction
%  		   -(1) - fwhm in focus (thus the resolution AND the angular point of the parabol) (sigma_0)
%          -(2) - depth of focus (d)
%		   -(3) - shifting of the parabol at the z-scale (arbitrary shift for computation)
%		   -(4) - gradient (m) of the straight line (mx+n) that biases the fwhm (no offset n) (m)
%		   -(5) - real vertex of the biased parabola
%		   -(6) - fwhm offset (usually negative)
%		   -(7) - lens shift
%	fwhmConfY - parameters for the fwhm curve in y direction (NOTE: 6 and 7 are 0 here!)
%	axialRange - defines the axial range of the data that is maximally necessary [lowerBound upperBound] in nm
%	minNumVertexFrames - minimal number of vertex frames that is required to start estimation of the fwhm curves
%	doOnlineCalibration - if 0 use predefined fwhm curve, if 1 do online calibration and estimate the curve
%	estimateFwhmOffset - if 0 fwhm offset must be provided and is fixed, if 1 offset will be estimated
%	doCombinedCLensFit - only important for evaluation (should be 1)
%	computeZPosition - 0 -> do not compute the axial position
%						1 -> already compute the axial position
%	returnCode: 
%					1 -> all cool
%					-1 -> didn't find enough vertex frames
%	control:	some verification cell array that contains certain measures to check if an image can provide correct results
%			1 - # of photons [mean std]
%			2 - peak to background ratio (note that the number of backgrund photons is used as per pixel such that even if the region size was changed the
%			values are comparable) [mean std]
%			3 - # focal place crossings [x y]
%			4 - lateral movement (euclidean dist) [d1 d2 d3 ..]
%			5 - fwhm at vertex [combined onlyForXCurve onlyForYCurve]
%			6 - track lengths (NOTE: it counts anything as a single track that is separated by a nan value) [l1 l2 l3 ...]
%			7 - estimated fwhmOffset [median mean]
%			8 - difference between the fwhmX and fwhmY at the vertices [medianX meanX medianY meanY]
%			9 - number of photons at vertex points [meanX stdX meanY stdY]
%			10 - peak to background ratio at vertex points [meanX stdX meanY stdY]
%			11 - cell array with some internal data that can be used to identify which tracks contributed to the calibration via their vertex points
%				(1) trFwhmX
%				(2) trFwhmY
%				(3) posX
%				(4) posY
%				(5) numPhotons
%				(6) bgRatio
%				(7) frameIdx
%				(8) vertexIdxX
%				(9) vertexIdxY
%				(10) fwhmCurveCoord
%				(11) fwhmRefCurveCoord
function [p fwhmConfX fwhmConfY fwhmCurveCoord fwhmRefCurveCoord returnCode control] = mda2(im, dv, trackingParams, plots, EMGain, p, fwhmConfX, fwhmConfY, axialRange, minNumVertexFrames, doOnlineCalibration, estimateFwhmOffset, doCombinedCLensFit, computeZPosition)%, vertexDistance)

%% constructor
returnCode = 1;
control = {};

% a variable for plotting purpose such that the fwhm curves have the same fwhm at their vertex point
xyFwhmOffset = 0;%fwhmConfX(6)%-0.508025338915882;%-43;%0;

%% parameters/variables

% cell array holding gConfs for each frame
% parameters:
% -gConf	- vector of parameters defining the gaussian:
%      -(1) - center x position (in fitting image frame)
%      -(2) - center y position (in fitting image frame)
%      -(3) - fwhm x
%      -(4) - fwhm y
%      -(5) - height of the peak
%      -(6) - offset
%      -(7) - center x position (in full image)
%      -(8) - center y position	(in full image)
%	   -(9) - turning angle theta if an elliptical gaussian is used
%	   -(10) - relative z position (using difference function)
%	   -(11) - relative z position (using mean z value of the two fwhm curves separately)
if isempty(p)
	p = {};
end

% if no trackingParams are provided take this one
if isempty(trackingParams)
	fprintf('WARNING: no trackingParams were provided, switching to default ones\n');
	
	trackingParams{1}=3;% maxVoxelDisplacement
	trackingParams{2}=0.5;% synapseDistInterval
	trackingParams{3}=[200, 200; 3000 3000];% minMaxFWHM
	trackingParams{4}=0;% gType
	trackingParams{5}=0;%8;% peakThreshRel
	trackingParams{6}=25;% minNumPhotonsPerPeak
	trackingParams{7}=4;% fitQualityThresholdFactor
	trackingParams{8}=2/3;% gaussianPart2Fit
	trackingParams{9}=0.1;% intensityOutliers
	trackingParams{10}=[0 1 0 0 1];% reasonabilityChecks
	trackingParams{11}=3;% subImSizefactor
	trackingParams{12}=1;% useWatershedTransformRegions
	trackingParams{13}=1;% use weighted LSQ or not
	trackingParams{14}=[400 400];% initFwhm
	trackingParams{15}=1.5;% watershedRegionControlAxisSizeFactor
	trackingParams{16}=1;% gaussianOffsetConstraints, 1: std-bounded, 2: fixed or 3: unbounded gaussian offset
	trackingParams{17}=0;% removeAllParticlesTooCloseToOthers 0: no 1: yes
	% setting: useImTopHatFilter = [iDiffRel, dTop, dBrim] == [relative intensity difference (as a factor for the std), diameter of top and brim]
	trackingParams{18}=[0.5 7 11];%[0.5 9 13];%[];% useImTopHatFilter
end

% acquisition parameters
% EMGain = 64;
CCDSensitivity = 12.4; % using Andor iXionEM-897 CCD camera

fwhmCurveCoord = [];
fwhmRefCurveCoord = [];

%% main

%%% find all possible occurences of particles in the 2D image slices

% if we want to compute all from scratch
if isempty(p)
	
	fprintf('find all possible occurences of particles in the 2D image slices..\n');
	
	% loop over frames
	for fNum = 1:size(im, 3)

		% TODO: check if the fitting result can be increased (and thus the axial position) because if i plot the results it looks like the watershed algorithm does
		% not leave much pixels for the fitting
		p{fNum} = locateParticles2DBetha(im(:,:, fNum), dv, trackingParams, [], 0, CCDSensitivity, EMGain, [], []);

		if size(p{fNum}, 1) == 0
			fprintf('EMPTY gConfs!!!!!!\n');
		end

		% so far we just want to track a single particle
		if size(p{fNum}, 1) > 1
			% take only first one and warn
			p{fNum} = p{fNum}(1,:);
			fprintf('\nWARNING: frame %d contains more than one particle - taking first one!\n\n', fNum);
			
		end

		% if we want to plot the fitting results over the image
		if plots(1)

			% compute mask covering the fitted qds
			IMask = zeros(size(im, 1), size(im, 2));

			% loop over particles and set the according number at pixels of the same region in the mask
			for pNum = 1 : size(p{fNum}, 1)

				% compute gaussian
				if trackingParams{4} == 2
					null = gaussian2DEllipse([p{fNum}(pNum, 7), p{fNum}(pNum, 8), p{fNum}(pNum, 3), p{fNum}(pNum, 4), 1, 0, 0, 0, p{fNum}(pNum, 9)], dv, size(im, 2), size(im, 1));
				else
					null = gaussian2D([p{fNum}(pNum, 7), p{fNum}(pNum, 8), p{fNum}(pNum, 3), p{fNum}(pNum, 4), 1, 0], dv, size(im, 2), size(im, 1));
				end

				% cut the gaussian at a certain sigma-interval and insert the pNum as the index for the
				% region
				% NOTE: this computation is actually independent of sigma
				null(null<exp(-0.5)) = 0;
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
			overlay(imresize(im(:,:,fNum), 1.5, 'nearest'), imresize(IMaskBinary, 1.5, 'nearest'), [0 1 0], 0.2, sprintf('#%d', fNum));
		end

	end
end

%%% do online calibration of the axial transformation function

% always store reference calibration curves
fwhmRefConfX = fwhmConfX;
fwhmRefConfY = fwhmConfY;

% if online calibration shall be performed
if doOnlineCalibration
% % if there are no calibration curves provided
% if isnan(fwhmConfX(1)) || isnan(fwhmConfY(1))
	
	fprintf('\n*****************************\n* Doing online calibration! *\n*****************************\n\n');

	trFwhmX = [];
	trFwhmY = [];
	numPhotons = [];
	bgRatio = [];
	posX = [];
	posY = [];
	
	% get number of pixels in the image
	% NOTE: by now it is not checked if this number is really correct for each particle because if a particle is close to the border this value might be
	% reduced, however using our manual regions in mdaGui right now this is not critical as they are usually centered
	if mod(trackingParams{11}, 2) == 0
		numPixels = (trackingParams{11}-1)^2;
	else
		numPixels = (trackingParams{11})^2;
	end
	
	%% check if we get just particle streams (first version) or if there are also the indexes provided (as required by mdaGui)
	if size(p, 2) == 1
		indexed = 0;
	else
		indexed = 1;
		frameIdx = [];
	end
		
	% loop over frames and get the fwhm for the trace
	for fNum = 1:size(p, 1)%numel(p)%size(im, 3)

		% if there was no particle
		if size(p{fNum}, 1) == 0

			if ~indexed
				% set fwhm to nan
				trFwhmX(fNum) = nan;
				trFwhmY(fNum) = nan;
				posX(fNum) = nan;
				posY(fNum) = nan;
				numPhotons(fNum) = nan;
				bgRatio(fNum) = nan;
			else
				% add nan value
				trFwhmX = [trFwhmX nan];
				trFwhmY = [trFwhmY nan];
				posX = [posX nan];
				posY = [posY nan];
				numPhotons = [numPhotons nan];
				bgRatio = [bgRatio nan];
				frameIdx = [frameIdx nan];
			end

			% and skip
			continue;
		end

		% loop over particles
		for pNum = 1:size(p{fNum}, 1)

			% so far we just want to track a single particle
			if pNum > 1
				% thus warn
				fprintf('\nWARNING: frame %d contains more than one particle - ignoring!\n\n', fNum);
				% and skip
				continue;
			end

			% collect fwhms for one particle
			if ~indexed
				trFwhmX(fNum) = p{fNum}(pNum, 3)*dv(1);
				trFwhmY(fNum) = p{fNum}(pNum, 4)*dv(2);
				posX(fNum) = p{fNum}(pNum, 7);
				posY(fNum) = p{fNum}(pNum, 8);
				numPhotons(fNum) = p{fNum}(pNum, 5);
				bgRatio(fNum) = p{fNum}(pNum, 5)/(p{fNum}(pNum, 6)/numPixels);
			else
				if isempty(trFwhmX) || isempty(p{fNum-1, 1})
					% just add the value because its either the first or after a nan value that was added above because there were no particles in the frame before
					trFwhmX = [trFwhmX p{fNum, 1}(pNum, 3)*dv(1)];
					trFwhmY = [trFwhmY p{fNum, 1}(pNum, 4)*dv(2)];
					posX = [posX p{fNum, 1}(pNum, 7)];
					posY = [posY p{fNum, 1}(pNum, 8)];
					numPhotons = [numPhotons p{fNum}(pNum, 5)];
					bgRatio = [bgRatio p{fNum}(pNum, 5)/(p{fNum}(pNum, 6)/numPixels)];
					frameIdx = [frameIdx p{fNum, 2}];
				elseif p{fNum, 2}-p{fNum-1, 2}>1
					% add nan AND the current value
					trFwhmX = [trFwhmX nan p{fNum, 1}(pNum, 3)*dv(1)];
					trFwhmY = [trFwhmY nan p{fNum, 1}(pNum, 4)*dv(2)];
					posX = [posX nan p{fNum, 1}(pNum, 7)];
					posY = [posY nan p{fNum, 1}(pNum, 8)];
					numPhotons = [numPhotons nan p{fNum}(pNum, 5)];
					bgRatio = [bgRatio nan p{fNum}(pNum, 5)/(p{fNum}(pNum, 6)/numPixels)];
					frameIdx = [frameIdx nan p{fNum, 2}];
				elseif p{fNum, 2}-p{fNum-1, 2}==1
					% also just add because its the directly following particle
					trFwhmX = [trFwhmX p{fNum, 1}(pNum, 3)*dv(1)];
					trFwhmY = [trFwhmY p{fNum, 1}(pNum, 4)*dv(2)];
					posX = [posX p{fNum, 1}(pNum, 7)];
					posY = [posY p{fNum, 1}(pNum, 8)];
					numPhotons = [numPhotons p{fNum}(pNum, 5)];
					bgRatio = [bgRatio p{fNum}(pNum, 5)/(p{fNum}(pNum, 6)/numPixels)];
					frameIdx = [frameIdx p{fNum, 2}];
				else
					fprintf('\nWARNING: frame %d could not be classified - exiting!\n\n', fNum);
					clear fwhmConfX;
					p{fNum, 2}-p{fNum-1, 2}
					p{fNum, 2}
					p{fNum-1, 2}
					return
				end
			end

% 			numPhotons = [numPhotons p{fNum}(pNum, 5)];
% 			bgRatio = [bgRatio p{fNum}(pNum, 5)/(p{fNum}(pNum, 6)/numPixels)];

% 			% if there is no particle
% 			if pNum == 0
% 				% set fwhm to nan
% 				trFwhmX(fNum) = nan;
% 				trFwhmY(fNum) = nan;
% 			else
% 				% collect fwhms for one particle
% 				trFwhmX(fNum) = p{fNum}(pNum, 3)*dv(1);
% 				trFwhmY(fNum) = p{fNum}(pNum, 4)*dv(2);
% 			end
		end
	end

	% save number of photons
	control{1} = [nanmean(numPhotons) nanstd(numPhotons)];
	% save peak to background ratio
	control{2} = [nanmean(bgRatio) nanstd(bgRatio)];
	% save track lengths (NOTE: it counts anything as a single track that is separated by a nan value)
	null = diff(find(isnan([nan trFwhmX nan])))-1;
	null(null==0) = [];
	control{6} = null;
	
	%% estimate lateral movement
	
	latMov = [];
	posX2 = [nan posX nan];
	posY2 = [nan posY nan];
	% find start index of tracks
	trIndexes = find(isnan([nan posX nan]));

	% loop over tracks
	for i = 1:length(trIndexes)-1
		% get single track
		singleTrackX = posX2(trIndexes(i)+1:trIndexes(i+1)-1);
		singleTrackY = posY2(trIndexes(i)+1:trIndexes(i+1)-1);
		
		if isempty(singleTrackX)
			continue;
		end
		
		% compute lateral movement
		latMov = [latMov sqrt((max(singleTrackX)-min(singleTrackX))^2 + (max(singleTrackY)-min(singleTrackY))^2)];
	end
	
	% save lateral movement
	control{4} = latMov;
	
%	% loop over frames and get the fwhm for the trace
% 	for fNum = 1:numel(p)%size(im, 3)
% 
% 		% if there was no particle
% 		if size(p{fNum}, 1) == 0
% 
% 			% set fwhm to nan
% 			trFwhmX(fNum) = nan;
% 			trFwhmY(fNum) = nan;
% 
% 			% and skip
% 			continue;
% 		end
% 
% 		% loop over particles
% 		for pNum = 1:size(p{fNum}, 1)
% 
% 			% so far we just want to track a single particle
% 			if pNum > 1
% 				% thus warn
% 				fprintf('\nWARNING: frame %d contains more than one particle - ignoring!\n\n', fNum);
% 				% and skip
% 				continue;
% 			end
% 
% 			% collect fwhms for one particle
% 			trFwhmX(fNum) = p{fNum}(pNum, 3)*dv(1);
% 			trFwhmY(fNum) = p{fNum}(pNum, 4)*dv(2);
% 
% % 			% if there is no particle
% % 			if pNum == 0
% % 				% set fwhm to nan
% % 				trFwhmX(fNum) = nan;
% % 				trFwhmY(fNum) = nan;
% % 			else
% % 				% collect fwhms for one particle
% % 				trFwhmX(fNum) = p{fNum}(pNum, 3)*dv(1);
% % 				trFwhmY(fNum) = p{fNum}(pNum, 4)*dv(2);
% % 			end
% 		end
% 	end

	
% 	figure;
% 	hold on;
% 	title(sprintf('plot fwhm x against fwhm y'));
% 	plot(trFwhmX, trFwhmY, 'b.', 'LineStyle', 'none');
% 	xlim([200 600])
% 	ylim([200 600])
% 	hold off;
% 
% 	figure;
% 	hold on;
% 	title('fwhm course and vertex detection');
% 	plot(1:length(trFwhmX), trFwhmX+xyFwhmOffset, 'b.');
% 	plot(1:length(trFwhmY), trFwhmY, 'r.');
% % 	hold off;
	
	% collect possible vertex points using structural, logical properties of the two curves:
	%	-logical threshold saying that a vertex point must be surrounded by two higher fwhm points and the same window in the other fwhm curve must have constantly
	%	decreasing or increasing values
	%	-quantil threshold
	fprintf('\nTODO: there must be further constraints because we usually underestimate the vertex value because lower outliers are too easily included!\n\n');
	vertexIdxX = [];
	vertexIdxY = [];
	
	%%% THREE NEIGHBORHOOD!!!!!!!!!
	% loop over trajectory and scan for vertex points
	for fNum = 2:length(trFwhmX)-1
	
		% check for fwhm vertex points
		% a vertex point must be surrounded by two higher fwhm points and the same window in the other fwhm curve must have constantly decreasing or increasing values
		
% 		-alternativ die robustheit steigern indem die 5 punkt umgebung eines vertex ausgenutzt wird
% 		-nur die besten 20% der vertexe als welche nehmen von denen dann l und r berechnet werden
		
		
		% xfwhm
		if trFwhmX(fNum-1) > trFwhmX(fNum) && trFwhmX(fNum+1) > trFwhmX(fNum)...
				&& ((trFwhmY(fNum-1) < trFwhmY(fNum) && trFwhmY(fNum) < trFwhmY(fNum+1))...
				|| (trFwhmY(fNum-1) > trFwhmY(fNum) && trFwhmY(fNum) > trFwhmY(fNum+1)))...
				&& ~(isnan(trFwhmX(fNum-1)) || isnan(trFwhmX(fNum)) || isnan(trFwhmX(fNum+1)) || isnan(trFwhmY(fNum-1)) || isnan(trFwhmY(fNum)) || isnan(trFwhmY(fNum+1)))
			
			% if the fwhm offset is known
			if ~estimateFwhmOffset
				% if the point has a lower value than the value in the other fwhm curve
				if trFwhmX(fNum)+fwhmConfX(6) < trFwhmY(fNum)
					% finally store the point as a vertex point
					vertexIdxX = [vertexIdxX fNum];
				end
			else
				% always store the point as a vertex point
				vertexIdxX = [vertexIdxX fNum];
			end
		end

		% yfwhm
		if trFwhmY(fNum-1) > trFwhmY(fNum) && trFwhmY(fNum+1) > trFwhmY(fNum)...
				&& ((trFwhmX(fNum-1) < trFwhmX(fNum) && trFwhmX(fNum) < trFwhmX(fNum+1))...
				|| (trFwhmX(fNum-1) > trFwhmX(fNum) && trFwhmX(fNum) > trFwhmX(fNum+1)))...
				&& ~(isnan(trFwhmX(fNum-1)) || isnan(trFwhmX(fNum)) || isnan(trFwhmX(fNum+1)) || isnan(trFwhmY(fNum-1)) || isnan(trFwhmY(fNum)) || isnan(trFwhmY(fNum+1)))

			% if the fwhm offset is known
			if ~estimateFwhmOffset
				% if the point has a lower value than the value in the other fwhm curve
				if trFwhmY(fNum) < trFwhmX(fNum)+fwhmConfX(6)
					% finally store the point as a vertex point
					vertexIdxY = [vertexIdxY fNum];
				end
			else
				% always store the point as a vertex point
				vertexIdxY = [vertexIdxY fNum];
			end
		end

	end

	%%% FIVE NEIGHBORHOOD!!!!!!!!!
% 	% loop over trajectory and scan for vertex points
% 	% NOTE: here the 5 neighborhood is analyzed instead of just three
% 	for fNum = 3:length(trFwhmX)-2
% 	
% 		% check for fwhm vertex points
% 		% a vertex point must be surrounded by four higher fwhm points and the same window in the other fwhm curve must have constantly decreasing or increasing values
% 		
% % 		% xfwhm
% % 		if trFwhmX(fNum-2) > trFwhmX(fNum-1) && trFwhmX(fNum-1) > trFwhmX(fNum) && trFwhmX(fNum+1) > trFwhmX(fNum) && trFwhmX(fNum+2) > trFwhmX(fNum+1)...
% % 				&& ((trFwhmY(fNum-2) < trFwhmY(fNum-1) && trFwhmY(fNum-1) < trFwhmY(fNum) && trFwhmY(fNum) < trFwhmY(fNum+1) && trFwhmY(fNum+1) < trFwhmY(fNum+2))...
% % 				|| (trFwhmY(fNum-2) > trFwhmY(fNum-1) && trFwhmY(fNum-1) > trFwhmY(fNum) && trFwhmY(fNum) > trFwhmY(fNum+1) && trFwhmY(fNum+1) > trFwhmY(fNum+2)))...
% % 				&& ~(isnan(trFwhmX(fNum-2)) || isnan(trFwhmX(fNum-1)) || isnan(trFwhmX(fNum)) || isnan(trFwhmX(fNum+1)) || isnan(trFwhmX(fNum+2))...
% % 				|| isnan(trFwhmY(fNum-2)) || isnan(trFwhmY(fNum-1)) || isnan(trFwhmY(fNum)) || isnan(trFwhmY(fNum+1)) || isnan(trFwhmY(fNum+2)))
% 
% 		% xfwhm
% 		% NOTE: here the linear change is not enforced
% 		if trFwhmX(fNum-2) > trFwhmX(fNum) && trFwhmX(fNum-1) > trFwhmX(fNum) && trFwhmX(fNum+1) > trFwhmX(fNum) && trFwhmX(fNum+2) > trFwhmX(fNum)...
% 				&& ((trFwhmY(fNum-2) < trFwhmY(fNum) && trFwhmY(fNum-1) < trFwhmY(fNum) && trFwhmY(fNum) < trFwhmY(fNum+1) && trFwhmY(fNum) < trFwhmY(fNum+2))...
% 				|| (trFwhmY(fNum-2) > trFwhmY(fNum) && trFwhmY(fNum-1) > trFwhmY(fNum) && trFwhmY(fNum) > trFwhmY(fNum+1) && trFwhmY(fNum) > trFwhmY(fNum+2)))...
% 				&& ~(isnan(trFwhmX(fNum-2)) || isnan(trFwhmX(fNum-1)) || isnan(trFwhmX(fNum)) || isnan(trFwhmX(fNum+1)) || isnan(trFwhmX(fNum+2))...
% 				|| isnan(trFwhmY(fNum-2)) || isnan(trFwhmY(fNum-1)) || isnan(trFwhmY(fNum)) || isnan(trFwhmY(fNum+1)) || isnan(trFwhmY(fNum+2)))
% 			
% 			% if the fwhm offset is known
% 			if ~estimateFwhmOffset
% 				% if the point has a lower value than the value in the other fwhm curve
% 				if trFwhmX(fNum)+fwhmConfX(6) < trFwhmY(fNum)
% 					% finally store the point as a vertex point
% 					vertexIdxX = [vertexIdxX fNum];
% 				end
% 			else
% 				% always store the point as a vertex point
% 				vertexIdxX = [vertexIdxX fNum];
% 			end
% 		end
% 
% % 		% yfwhm
% % 		if trFwhmY(fNum-2) > trFwhmY(fNum-1) && trFwhmY(fNum-1) > trFwhmY(fNum) && trFwhmY(fNum+1) > trFwhmY(fNum) && trFwhmY(fNum+2) > trFwhmY(fNum+1)...
% % 				&& ((trFwhmX(fNum-2) < trFwhmX(fNum-1) && trFwhmX(fNum-1) < trFwhmX(fNum) && trFwhmX(fNum) < trFwhmX(fNum+1) && trFwhmX(fNum+1) < trFwhmX(fNum+2))...
% % 				|| (trFwhmX(fNum-2) > trFwhmX(fNum-1) && trFwhmX(fNum-1) > trFwhmX(fNum) && trFwhmX(fNum) > trFwhmX(fNum+1) && trFwhmX(fNum+1) > trFwhmX(fNum+2)))...
% % 				&& ~(isnan(trFwhmX(fNum-2)) || isnan(trFwhmX(fNum-1)) || isnan(trFwhmX(fNum)) || isnan(trFwhmX(fNum+1)) || isnan(trFwhmX(fNum+2))...
% % 				|| isnan(trFwhmY(fNum-2)) || isnan(trFwhmY(fNum-1)) || isnan(trFwhmY(fNum)) || isnan(trFwhmY(fNum+1)) || isnan(trFwhmY(fNum+2)))
% 
% 		% yfwhm
% 		% NOTE: here the linear change is not enforced
% 		if trFwhmY(fNum-2) > trFwhmY(fNum) && trFwhmY(fNum-1) > trFwhmY(fNum) && trFwhmY(fNum+1) > trFwhmY(fNum) && trFwhmY(fNum+2) > trFwhmY(fNum)...
% 				&& ((trFwhmX(fNum-2) < trFwhmX(fNum) && trFwhmX(fNum-1) < trFwhmX(fNum) && trFwhmX(fNum) < trFwhmX(fNum+1) && trFwhmX(fNum) < trFwhmX(fNum+2))...
% 				|| (trFwhmX(fNum-2) > trFwhmX(fNum) && trFwhmX(fNum-1) > trFwhmX(fNum) && trFwhmX(fNum) > trFwhmX(fNum+1) && trFwhmX(fNum) > trFwhmX(fNum+2)))...
% 				&& ~(isnan(trFwhmX(fNum-2)) || isnan(trFwhmX(fNum-1)) || isnan(trFwhmX(fNum)) || isnan(trFwhmX(fNum+1)) || isnan(trFwhmX(fNum+2))...
% 				|| isnan(trFwhmY(fNum-2)) || isnan(trFwhmY(fNum-1)) || isnan(trFwhmY(fNum)) || isnan(trFwhmY(fNum+1)) || isnan(trFwhmY(fNum+2)))
% 
% 			% if the fwhm offset is known
% 			if ~estimateFwhmOffset
% 				% if the point has a lower value than the value in the other fwhm curve
% 				if trFwhmY(fNum) < trFwhmX(fNum)+fwhmConfX(6)
% 					% finally store the point as a vertex point
% 					vertexIdxY = [vertexIdxY fNum];
% 				end
% 			else
% 				% always store the point as a vertex point
% 				vertexIdxY = [vertexIdxY fNum];
% 			end
% 		end
% 
% 	end

% 	% estimate vertex point of the fwhm curves separately
% 	% if the fwhm offset is known
% 	if ~estimateFwhmOffset
% 		vertexFwhmX = median(trFwhmX(vertexIdxX)+fwhmConfX(6));
% 		vertexFwhmXMean = mean(trFwhmX(vertexIdxX)+fwhmConfX(6));
% 	% if the fwhm offset needs to be estimated
% 	else
% 		vertexFwhmX = median(trFwhmX(vertexIdxX));
% 		vertexFwhmXMean = mean(trFwhmX(vertexIdxX));
% 	end
% 	vertexFwhmY = median(trFwhmY(vertexIdxY));
% 	vertexFwhmYMean = mean(trFwhmY(vertexIdxY));

	numVertexPointsX  = length(vertexIdxX)
	numVertexPointsY  = length(vertexIdxY)
	
	% estimate the true vertex point
	% if the fwhm offset is known
	if ~estimateFwhmOffset
		vertexFwhm = median([trFwhmX(vertexIdxX)+fwhmConfX(6) trFwhmY(vertexIdxY)]);
% 		vertexFwhm = mean([trFwhmX(vertexIdxX)+fwhmConfX(6) trFwhmY(vertexIdxY)]);
	% if the fwhm offset needs to be estimated
	else
		% estimate fwhm offset and store it in the fwhmConfX
		fwhmConfX(6) = -(median(trFwhmX(vertexIdxX))-median(trFwhmY(vertexIdxY)));
% 		fwhmConfX(6) = -mean(trFwhmX(vertexIdxX)-trFwhmY(vertexIdxY));
		
		vertexFwhm = median([trFwhmX(vertexIdxX)+fwhmConfX(6) trFwhmY(vertexIdxY)]);
%		vertexFwhm = mean([trFwhmX(vertexIdxX)+fwhmConfX(6) trFwhmY(vertexIdxY)]);
	end
	
% 	% plot vertex points
% 	plot(vertexIdxX, trFwhmX(vertexIdxX)+xyFwhmOffset, 'g+');
% 	plot(vertexIdxY, trFwhmY(vertexIdxY), 'k+');
% 	hold off;
% 
% 	% plot the difference in the fwhm between the vertex points and the fwhm of the orthogonal direction
% 	figure;
% 	hold on;
% 	title('fwhm difference at vertex points');
% 	plot(trFwhmY(vertexIdxX)-(trFwhmX(vertexIdxX)+fwhmConfX(6)), 'b.', 'LineStyle', 'none');
% 	plot((trFwhmX(vertexIdxY)+fwhmConfX(6))-trFwhmY(vertexIdxY), 'r.', 'LineStyle', 'none');
% 	hold off;
% 	
% 	figure;
% 	hold on;
% 	title('plot VERTEX fwhm x against fwhm y');
% 	plot(trFwhmX(vertexIdxX)+fwhmConfX(6), trFwhmY(vertexIdxX), 'b.', 'LineStyle', 'none');
% 	plot(trFwhmX(vertexIdxY)+fwhmConfX(6), trFwhmY(vertexIdxY), 'r.', 'LineStyle', 'none');
% 	xlim([200 600])
% 	ylim([200 600])
% 	hold off;
% 	
% 	debugStop = 1;

% 	% estimate difference to the other fwhm curve at the vertex point of each fwhm curve
% 	vertexFwhmLeft = median(trFwhmY(vertexIdxX)-(trFwhmX(vertexIdxX)+fwhmConfX(6)))
% % 	vertexFwhmLeft = mean(trFwhmY(vertexIdxX)-(trFwhmX(vertexIdxX)+fwhmConfX(6)));
% 	vertexFwhmRight = median((trFwhmX(vertexIdxY)+fwhmConfX(6))-trFwhmY(vertexIdxY))
% % 	vertexFwhmRight = mean((trFwhmX(vertexIdxY)+fwhmConfX(6))-trFwhmY(vertexIdxY));

	%% restrict the vertex points to those that are close to the mean value such that only "good" vertices are used for the computation of the coordinates left
	%% and right of the vertex
	
% 	fprintf('ATTENTION: performing restriction of vertex points!\n');
% 	% NOTE: this restriction is reasonable as the mean estimation for the vertex is quite reliable
% 	vertexIdxXRestricted = [];
% 	vertexIdxYRestricted = [];
% 	% compute acceptable deviation
% % % 	stdVertices = mad([trFwhmX(vertexIdxX)+fwhmConfX(6) trFwhmY(vertexIdxY)], 1)
% % 	stdVertices = std([trFwhmX(vertexIdxX)+fwhmConfX(6) trFwhmY(vertexIdxY)])
% 	if numel(vertexIdxX)>0
% 		stdVerticesX = std(trFwhmX(vertexIdxX)+fwhmConfX(6));
% 	else
% 		stdVerticesX = 0;
% 	end
% 	if numel(vertexIdxY)>0
% 		stdVerticesY = std(trFwhmY(vertexIdxY));
% 	else
% 		stdVerticesY = 0;
% 	end
% 	
% 	for vNum = 1:length(vertexIdxX)
% 		
% 		% if the current vertex estimate deviates too much (right now much more as the median absolute deviation) from the mean than it is rejected for the new list
% 		if abs(vertexFwhm-(trFwhmX(vertexIdxX(vNum))+fwhmConfX(6))) < stdVerticesX
% 			% store again in new list
% 			vertexIdxXRestricted = [vertexIdxXRestricted vertexIdxX(vNum)];
% 		end
% 	end
% 	for vNum = 1:length(vertexIdxY)
% 		
% 		% if the current vertex estimate deviates too much (right now much more as the median absolute deviation) from the mean than it is rejected for the new list
% 		if abs(vertexFwhm-(trFwhmY(vertexIdxY(vNum)))) < stdVerticesY
% 			% store again in new list
% 			vertexIdxYRestricted = [vertexIdxYRestricted vertexIdxY(vNum)];
% 		end
% 	end
% 
% 	numVertexPointsXRestricted = length(vertexIdxXRestricted)
% 	numVertexPointsYRestricted = length(vertexIdxYRestricted)
% 	
% 	% recompute vertex Fwhm
% 	vertexFwhm = median([trFwhmX(vertexIdxX)+fwhmConfX(6) trFwhmY(vertexIdxY)]);
% 	% overwrite vertex indices
% 	vertexIdxX = vertexIdxXRestricted;
% 	vertexIdxY = vertexIdxYRestricted;

	% store number of vertex points
	control{3} = [length(vertexIdxX) length(vertexIdxY)];
	% store fwhm at vertex [combined onlyForXCurve onlyForYCurve]
	control{5} = [vertexFwhm median(trFwhmX(vertexIdxX)+fwhmConfX(6)) median(trFwhmY(vertexIdxY))];
	% store the difference between the fwhmX and fwhmY at the vertices [medianX meanX medianY meanY]
	control{8} = [median(trFwhmY(vertexIdxX)-(trFwhmX(vertexIdxX)+fwhmConfX(6))) mean(trFwhmY(vertexIdxX)-(trFwhmX(vertexIdxX)+fwhmConfX(6))) median((trFwhmX(vertexIdxY)+fwhmConfX(6))-trFwhmY(vertexIdxY)) mean((trFwhmX(vertexIdxY)+fwhmConfX(6))-trFwhmY(vertexIdxY))];
	% store the estimated fwhmOffset [median mean]
	control{7} = [median(trFwhmX(vertexIdxX))-median(trFwhmY(vertexIdxY)) mean(trFwhmX(vertexIdxX))-mean(trFwhmY(vertexIdxY))];
	% save number of photons at vertex points [meanX stdX meanY stdY]
	control{9} = [nanmean(numPhotons(vertexIdxX)) nanstd(numPhotons(vertexIdxX)) nanmean(numPhotons(vertexIdxY)) nanstd(numPhotons(vertexIdxY))];
	% save peak to background ratio at vertex points [meanX stdX meanY stdY]
	control{10} = [nanmean(bgRatio(vertexIdxX)) nanstd(bgRatio(vertexIdxX)) nanmean(bgRatio(vertexIdxY)) nanstd(bgRatio(vertexIdxY))];
	
	%% skip curve estimation in case there are not enough vertex points found
	if length(vertexIdxX) < minNumVertexFrames || length(vertexIdxY) < minNumVertexFrames
		fprintf('WARNING: found unsufficient number of vertex frames (<%d), skipping curve estimation!\n', minNumVertexFrames);
		returnCode = -1;
		return;
	end
	
	%% construct the threee candidates
	% due to the two fwhm curves having the same curvature but being shifted at the axial axis we now already know three points of each curve which are:
	%	-the vertex point with the estimated fwhm (y) and zero axial (x) value
	%	-the points +-vertexDistance (x) and their according estimated fwhm (y) value, which, due to the similarity of the curves, are computed from the
	%	values in vertexDistance in the direction of the vertex pint of the other fwhm curve
	fwhmCurveCoord(1, 1:2) = [0 vertexFwhm]; % vertex point
% % 	fwhmCurveCoord(2, 1:2) = [-fwhmConfX(7) vertexFwhmLeft+vertexFwhm]; % estimated fwhm value in vertexDistance left of the vertex point
% % 	fwhmCurveCoord(3, 1:2) = [fwhmConfX(7) vertexFwhmRight+vertexFwhm]; % estimated fwhm value in vertexDistance right of the vertex point
% % 	fwhmCurveCoord(2, 1:2) = [-fwhmConfX(7) mean(trFwhmY(vertexIdxX))]; % estimated fwhm value in vertexDistance left of the vertex point
% % 	fwhmCurveCoord(3, 1:2) = [fwhmConfX(7) mean((trFwhmX(vertexIdxY)+fwhmConfX(6)))]; % estimated fwhm value in vertexDistance right of the vertex point
	fwhmCurveCoord(2, 1:2) = [-fwhmConfX(7) median(trFwhmY(vertexIdxX))]; % estimated fwhm value in vertexDistance left of the vertex point
	fwhmCurveCoord(3, 1:2) = [fwhmConfX(7) median((trFwhmX(vertexIdxY)+fwhmConfX(6)))]; % estimated fwhm value in vertexDistance right of the vertex point
% 	fwhmCurveCoord(2, 1:2) = [-fwhmConfX(7) median(trFwhmY(vertexIdxXRestricted))]; % estimated fwhm value in vertexDistance left of the vertex point
% 	fwhmCurveCoord(3, 1:2) = [fwhmConfX(7) median((trFwhmX(vertexIdxYRestricted)+fwhmConfX(6)))]; % estimated fwhm value in vertexDistance right of the vertex point

	fwhmCurveCoord

	if doCombinedCLensFit
		% compare to reference points at positions that should be the same in both curves due to the assumption that they are similar.
		% however, most probably one should not assume that for the reference curve but only as an assumption for our model that is applied later on because
		% else we add the deviation of the assumption twice: in the reference and the model
		fwhmRefCurveCoord(1, 1:2) = [0 biasedFwhm1DDiscrete(fwhmRefConfY, 0, dv(3))]; % vertex point
		fwhmRefCurveCoord(2, 1:2) = [-fwhmRefConfX(7) biasedFwhm1DDiscrete(fwhmRefConfY, -fwhmRefConfX(7), dv(3))]; % estimated fwhm value in vertexDistance left of the vertex point
		fwhmRefCurveCoord(3, 1:2) = [fwhmRefConfX(7) biasedFwhm1DDiscrete(fwhmRefConfY, fwhmRefConfX(7), dv(3))]; % estimated fwhm value in vertexDistance right of the vertex point
	else
		% compare to the reference points that at their real position from where we estimated it
		fwhmRefCurveCoord(1, 1:2) = [0 biasedFwhm1DDiscrete(fwhmRefConfY, 0, dv(3))]; % vertex point
		fwhmRefCurveCoord(2, 1:2) = [-fwhmRefConfX(7) biasedFwhm1DDiscrete(fwhmRefConfY, -fwhmRefConfX(7), dv(3))]; % estimated fwhm value in vertexDistance left of the vertex point
		% NOTE: use position zero here because the vertex of the y curve is shifted to the zero position
		fwhmRefCurveCoord(3, 1:2) = [fwhmRefConfX(7) biasedFwhm1DDiscrete(fwhmRefConfX, 0, dv(3))+fwhmConfX(6)]; % estimated fwhm value in vertexDistance right of the vertex point
	end
	fwhmRefCurveCoord

 	% store some internal data that can be used to identify which tracks contributed to the calibration via their vertex points
	control{11} = {trFwhmX; trFwhmY; posX; posY; numPhotons; bgRatio; frameIdx; vertexIdxX; vertexIdxY; fwhmCurveCoord; fwhmRefCurveCoord};
	
% 	% test with original y curve fit
% % 	fwhmCurveCoord(1, 1:2) = [70.734263729817087-62 275.4066]; % vertex point
% % 	fwhmCurveCoord(2, 1:2) = [70.734263729817087-62-18 314.5182]; % estimated fwhm value in vertexDistance left of the vertex point
% % 	fwhmCurveCoord(3, 1:2) = [70.734263729817087-62+18 405.0224]; % estimated fwhm value in vertexDistance right of the vertex point
% 	fwhmCurveCoord(1, 1:2) = [0 275.4066]; % vertex point
% 	fwhmCurveCoord(2, 1:2) = [-18 314.5182]; % estimated fwhm value in vertexDistance left of the vertex point
% 	fwhmCurveCoord(3, 1:2) = [18 405.0224]; % estimated fwhm value in vertexDistance right of the vertex point
% 	fwhmCurveCoord

	%%% solve the nonlinear equation system using constrained curve fitting

	% just for communication with the functions to be fitted
	global axialRes;
	axialRes = dv(3);
% 	global vertexPointCoord;
% 	vertexPointCoord = 22.9750;%fwhmCurveCoord(1, 1);%22.9750;%31.7092;
	
% 	init = [539 273 94-62 1.7];
% 	lowerBound = [300 100 0 0.1];
% 	upperBound = [1000 700 26*2 10];
% 	init = [539 273 1.7];
% 	init = [411 373 0.7];
% 	lowerBound = [300 100 0.1];
	init = [400 300 0.5];
	lowerBound = [1 1 -10];
	upperBound = [2000 2000 10];
	
	% fit the curve
	[var null residual] = lsqcurvefit(@biasedFwhm1DFunc, init, fwhmCurveCoord(1:3, 1)', fwhmCurveCoord(1:3, 2)', lowerBound, upperBound, optimset('Display', 'off'));
	var
residual
	% create estimated fwhm configuration
	% NOTE: the center of the particle is defined by the center of the (biased) fwhmY curve because this is the real particle center as the fwhmX is just shifted by the
	% lens. therefore, finally conf(5) defines the center!
	% NOTE: the fwhm offset is just added as parameter #6 and must be considered by the following analysis
% 	% NOTE: the two curves are centered arround their vertex distance and thus the x curve is already shifted according to the lens shift
% 	fwhmEstConfX = [var(1) var(2) -(fwhmConfX(7)/dv(3))/2 var(3) 0	fwhmConfX(6)	fwhmConfX(7)];
% 	fwhmEstConfY = [var(1) var(2) +(fwhmConfX(7)/dv(3))/2 var(3) 0	0				0];
% 	fwhmEstConfX = [var(1) var(2) -fwhmConfX(7)/dv(3)	var(3) 0	fwhmConfX(6)	fwhmConfX(7)];
	fwhmEstConfX = [var(1) var(2) -fwhmConfX(7)	var(3) 0	fwhmConfX(6)	fwhmConfX(7)];
	fwhmEstConfY = [var(1) var(2) 0				var(3) 0	0				0];

	% update the real position of the minima of the biased parabola
	fwhmEstConfY(5) = sqrt( (fwhmEstConfY(4)^2*fwhmEstConfY(2)^4) / (fwhmEstConfY(1)^2- (fwhmEstConfY(4)^2*fwhmEstConfY(2)^2) ));
	% correct the direction
	if fwhmEstConfY(4)>0
		fwhmEstConfY(5) = -fwhmEstConfY(5);
	end
	% add the axial shift and scale back to integer grid
	fwhmEstConfY(5) = fwhmEstConfY(5)/dv(3)+fwhmEstConfY(3);
	fwhmEstConfX(5) = fwhmEstConfY(5)-fwhmConfX(7);

	% finally really centralize to the real (biased) center of fwhmY
	fwhmEstConfX([3 5]) = fwhmEstConfX([3 5])-fwhmEstConfY(5);
	fwhmEstConfY([3 5]) = fwhmEstConfY([3 5])-fwhmEstConfY(5);

% 	% UNCOMMENT if the center shall be inbetween the two curve vertices
% 	fprintf('\n\nWARNING: moving center inbetween the two curves!!!\n\n');
% 	fwhmEstConfX([3 5]) = fwhmEstConfX([3 5])+((fwhmRefConfY(5)-fwhmRefConfX(5))/2);
% 	fwhmEstConfY([3 5]) = fwhmEstConfY([3 5])+((fwhmRefConfY(5)-fwhmRefConfX(5))/2);
	
% 	%%% solve the nonlinear equation system using symbolic math toolbox
% 	
% 	% create symbolic variables
% 	syms s0 d m
% 
% % 	% create three equations from three points
% % 	shift = 0;%-(70.734263729817087-62);
% % 	equ1 = sym('s0*sqrt(1+((fwhmCurveCoordCenterX-shift)/d)^2)+m*(fwhmCurveCoordCenterX-shift) = fwhmCurveCoordCenterY');
% % 	equ1 = subs(equ1, '[fwhmCurveCoordCenterX fwhmCurveCoordCenterY shift]', [fwhmCurveCoord(1, 1)*dv(3) fwhmCurveCoord(1, 2) shift]);
% % 	equ2 = sym('s0*sqrt(1+((fwhmCurveCoordLeftX-shift)/d)^2)+m*(fwhmCurveCoordLeftX-shift) = fwhmCurveCoordLeftY');
% % 	equ2 = subs(equ2, '[fwhmCurveCoordLeftX fwhmCurveCoordLeftY shift]', [fwhmCurveCoord(2, 1)*dv(3) fwhmCurveCoord(2, 2) shift]);
% % 	equ3 = sym('s0*sqrt(1+((fwhmCurveCoordRightX-shift)/d)^2)+m*(fwhmCurveCoordRightX-shift) = fwhmCurveCoordRightY');
% % 	equ3 = subs(equ3, '[fwhmCurveCoordRightX fwhmCurveCoordRightY shift]', [fwhmCurveCoord(3, 1)*dv(3) fwhmCurveCoord(3, 2) shift]);
% % % 	% add constraint that fwhmCurveCoord(1,:) is the vertex point by adding the equation saying that the first derivation at that point must be zero
% % % 	equ1 = sym('(s0/(d^2*sqrt(1+(fwhmCurveCoordCenterX/d)^2)))*fwhmCurveCoordCenterX+m = 0');
% % % 	equ1 = subs(equ1, '[fwhmCurveCoordCenterX]', [fwhmCurveCoord(1, 1)*dv(3)]);
% 
% % 	% create three equations from three points substituting for the shift on z
% % 	equ1 = sym('s0*sqrt(1+((fwhmCurveCoordCenterX+sqrt((m^2*d^4)/(s0^2-m^2*d^2)))/d)^2)+m*(fwhmCurveCoordCenterX+sqrt((m^2*d^4)/(s0^2-m^2*d^2))) = fwhmCurveCoordCenterY');
% % 	equ1 = subs(equ1, '[fwhmCurveCoordCenterX fwhmCurveCoordCenterY]', [fwhmCurveCoord(1, 1)*dv(3) fwhmCurveCoord(1, 2)]);
% % 	equ2 = sym('s0*sqrt(1+((fwhmCurveCoordLeftX+sqrt((m^2*d^4)/(s0^2-m^2*d^2)))/d)^2)+m*(fwhmCurveCoordLeftX+sqrt((m^2*d^4)/(s0^2-m^2*d^2))) = fwhmCurveCoordLeftY');
% % 	equ2 = subs(equ2, '[fwhmCurveCoordLeftX fwhmCurveCoordLeftY]', [fwhmCurveCoord(2, 1)*dv(3) fwhmCurveCoord(2, 2)]);
% % 	equ3 = sym('s0*sqrt(1+((fwhmCurveCoordRightX+sqrt((m^2*d^4)/(s0^2-m^2*d^2)))/d)^2)+m*(fwhmCurveCoordRightX+sqrt((m^2*d^4)/(s0^2-m^2*d^2))) = fwhmCurveCoordRightY');
% % 	equ3 = subs(equ3, '[fwhmCurveCoordRightX fwhmCurveCoordRightY]', [fwhmCurveCoord(3, 1)*dv(3) fwhmCurveCoord(3, 2)]);
% % % 	% add constraint that fwhmCurveCoord(1,:) is the vertex point by adding the equation saying that the first derivation at that point must be zero
% % % 	equ1 = sym('(s0/(d^2*sqrt(1+((fwhmCurveCoordCenterX+sqrt((m^2*d^4)/(s0^2-m^2*d^2)))/d)^2)))*(fwhmCurveCoordCenterX+sqrt((m^2*d^4)/(s0^2-m^2*d^2)))+m = 0');
% % % 	equ1 = subs(equ1, '[fwhmCurveCoordCenterX]', [fwhmCurveCoord(1, 1)*dv(3)]);
% 
% 	% create three equations from three points (original)
% 	equ1 = sym('s0*sqrt(1+(fwhmCurveCoordCenterX/d)^2)+m*fwhmCurveCoordCenterX = fwhmCurveCoordCenterY');
% 	equ1 = subs(equ1, '[fwhmCurveCoordCenterX fwhmCurveCoordCenterY]', [fwhmCurveCoord(1, 1)*dv(3) fwhmCurveCoord(1, 2)]);
% 	equ2 = sym('s0*sqrt(1+(fwhmCurveCoordLeftX/d)^2)+m*fwhmCurveCoordLeftX = fwhmCurveCoordLeftY');
% 	equ2 = subs(equ2, '[fwhmCurveCoordLeftX fwhmCurveCoordLeftY]', [fwhmCurveCoord(2, 1)*dv(3) fwhmCurveCoord(2, 2)]);
% 	equ3 = sym('s0*sqrt(1+(fwhmCurveCoordRightX/d)^2)+m*fwhmCurveCoordRightX = fwhmCurveCoordRightY');
% 	equ3 = subs(equ3, '[fwhmCurveCoordRightX fwhmCurveCoordRightY]', [fwhmCurveCoord(3, 1)*dv(3) fwhmCurveCoord(3, 2)]);
% % 	% add constraint that fwhmCurveCoord(1,:) is the vertex point by adding the equation saying that the first derivation at that point must be zero
% % 	equ1 = sym('(s0/(d^2*sqrt(1+(fwhmCurveCoordCenterX/d)^2)))*fwhmCurveCoordCenterX+m = 0');
% % 	equ1 = subs(equ1, '[fwhmCurveCoordCenterX]', [fwhmCurveCoord(1, 1)*dv(3)]);
% 
% % % 	equ1 = sym('s0*sqrt(1+(fwhmCurveCoordCenterX/d)^2)+m*fwhmCurveCoordCenterX = fwhmCurveCoordCenterY');
% % % % 	equ1 = subs(equ1, '[fwhmCurveCoordCenterX fwhmCurveCoordCenterY]', [fwhmCurveCoord(1, 1)*dv(3) fwhmCurveCoord(1, 2)]);
% % % 	equ1 = diff(equ1, 'fwhmCurveCoordCenterX')
% % % % 	equ1 = sym('(s0/(d^2*sqrt(1+(fwhmCurveCoordCenterX/d)^2)))*fwhmCurveCoordCenterX+m = 0');
% % % 	equ1 = subs(equ1, '[fwhmCurveCoordCenterX]', [fwhmCurveCoord(1, 1)*dv(3)]);
% 
% 	% solve linear equation system
%  	symVar = solve(equ1, equ2, equ3);%, equ4);
% 
% %  	symVar = solve('s0*sqrt(1+(subs(fwhmCurveCoordCenterX)/d)^2)+m*subs(fwhmCurveCoordCenterX) = subs(fwhmCurveCoordCenterY)'...
% %  		, 's0*sqrt(1+(subs(fwhmCurveCoordLeftX)/d)^2)+m*subs(fwhmCurveCoordLeftX) = subs(fwhmCurveCoordLeftY)'...
% %  		, 's0*sqrt(1+(subs(fwhmCurveCoordRightX)/d)^2)+m*subs(fwhmCurveCoordRightX) = subs(fwhmCurveCoordRightY)')
% % 	symVar = solve('s0*sqrt(1+(0/d)^2)+m*0 = 273.8229'...
% %  		, 's0*sqrt(1+((-360)/d)^2)+m*(-360) = 309.0834'...
% %  		, 's0*sqrt(1+(360/d)^2)+m*360 = 331.6805');
% 	
% % 	subs(symVar.s0)
% % 	subs(symVar.d)
% % 	subs(symVar.m)
% % 	symVar.s0
% % 	symVar.d
% % 	symVar.m
% 	
% 	% do reasonability check
% 	if isequal(subs(symVar.s0(1)), subs(symVar.s0(2))) && isequal(subs(symVar.m(1)), subs(symVar.m(2))) && isequal(subs(symVar.d(1)), subs(-symVar.d(2)))
% 		
% 		% and create estimated fwhm configuration
% 		% NOTE: the two curves are centered arround their vertex distance (anyway the absolute value doesn't matter)
% 		fwhmEstConfX = [subs(symVar.s0(1)) subs(symVar.d(1)) -(fwhmConfX(7)/dv(3))/2 subs(symVar.m(1))];
% 		fwhmEstConfY = [subs(symVar.s0(1)) subs(symVar.d(1)) +(fwhmConfX(7)/dv(3))/2 subs(symVar.m(1))];
% 	else
% 		fprintf('\nERROR: online calibration failed for reasonability check - exiting!\n\n');
% 		p = {};
% 		return;
% 	end

% 	%%% solve the nonlinear equation system using numerical solution
% 	global fwhmCurveCoord;
% 
% 	options=optimset('Display','iter', 'MaxFunEvals', 100000*4, 'MaxIter', 100000);
% 	[var, sValue] = fsolve(@equationSystem, [539 273 1.7], options)
% 	
% 	% create estimated fwhm configuration
% 	% NOTE: the two curves are centered arround their vertex distance (anyway the absolute value doesn't matter)
% 	fwhmEstConfX = [var(1) var(2) -(fwhmConfX(7)/dv(3))/2 var(3)];
% 	fwhmEstConfY = [var(1) var(2) +(fwhmConfX(7)/dv(3))/2 var(3)];
		
% 	% plot real and estimated fwhm curves
% 	figure;
% 	hold on;
% 	plot([-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3), biasedFwhm1DDiscrete(fwhmConfX, [-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3), dv(3)), 'b--');
% 	plot([-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3), biasedFwhm1DDiscrete(fwhmConfY, [-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3), dv(3)), 'r--');
% % 	plot([-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3), biasedFwhm1DDiscrete(fwhmEstConfX, [-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3), dv(3)), 'b');
% % 	plot([-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3), biasedFwhm1DDiscrete(fwhmEstConfY, [-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3), dv(3)), 'r');
% 	plot(fwhmCurveCoord(1:3, 1), fwhmCurveCoord(1:3, 2), 'go');
% 
% 	plot([-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3), biasedFwhm1DDiscrete(fwhmEstConfX, [-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3)-(fwhmConfX(7)/dv(3))/2, dv(3)), 'b');
% 	plot([-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3), biasedFwhm1DDiscrete(fwhmEstConfY, [-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3)-(fwhmConfX(7)/dv(3))/2, dv(3)), 'r');
% % 	plot(fwhmCurveCoord(1:3, 1)+[0 -36 +36]', fwhmCurveCoord(1:3, 2), 'go');
% 
% % 	plot([-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3), biasedFwhm1DDiscrete(fwhmConfY, -([-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3)-(70.734263729817087-62))+(fwhmConfX(7)/dv(3))/2, dv(3)), 'r--');
% % 	plot([-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3), biasedFwhm1DDiscrete(fwhmEstConfX, [-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3)-(fwhmConfX(7)/dv(3))/2, dv(3)), 'b');
% % 	plot([-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3), biasedFwhm1DDiscrete(fwhmEstConfY, [-2*fwhmConfX(7):2*fwhmConfX(7)]./dv(3)-(fwhmConfX(7)/dv(3))/2, dv(3)), 'r');
% % 	plot(fwhmCurveCoord(1:3, 1), fwhmCurveCoord(1:3, 2), 'go');
% 	hold off;
	
% 	biasedVertexFwhmYRef = biasedFwhm1DDiscrete(fwhmRefConfY, fwhmRefConfY(5), dv(3))

	fwhmRefConfX
	fwhmRefConfY

	% just for debugging we are copying the configuration that late
	fwhmConfX = fwhmEstConfX
	fwhmConfY = fwhmEstConfY

% 	% !!!!!!!UNCOMMENT!!!!!!!!! if you want to fake calibration results
%
% 	% NOTE: SHOULD ONLY BE USED FOR SYNTHETIC PSF. these are calibration results from analyzeRIM taken at zero depth with separated thus real clens fitting under normal conditions:
% 	% 	-removeFirstSliceInStack = 1
% 	% 	-doCombinedCLensFit = 0
% 	% 	-intensityThreshold = 0.2 % 0.12
% 	% 	-alignmentAlgorithm = 1
% 	% these values are used to simulate mdaSCRIPT2() results if always only the zero depth calibration is taken as the calibration curve. thus, we show that the
% 	% error becomes big in this case
% 	fprintf('\n\nWARNING: resetting calibration results to general calibration at zero depth!\n\n');
% 	fwhmConfX = [470.4750  381.2324   69.9002    0.9014   49.5146  -62.6313   11.3391];
% 	fwhmConfY = [259.4040  298.9264   59.7490   -0.0640   60.8537         0         0];
% 	% align all configurations to the real (biased) center of fwhmY
% 	fwhmConfX([3 5]) = fwhmConfX([3 5])-fwhmConfY(5);
% 	fwhmConfY([3 5]) = fwhmConfY([3 5])-fwhmConfY(5);
% 	% erase fwhmOffset
% 	fwhmConfX(6) = 0;
% 	% because these values were taken from acquisitions that had a dv(3) = 20 we need to adjust that to our resolution
% 	fwhmConfX([3 5 6 7]) = fwhmConfX([3 5 6 7]).*20;
% 	fwhmConfY([3 5 6 7]) = fwhmConfY([3 5 6 7]).*20;
% 	% overwrite the first four parameters of fwhmConfXRef because we want a really symmetric and perfekt aberrated PSF
% 	fwhmConfX([1 2 4]) = fwhmConfY([1 2 4]);
% 	fwhmConfX([3 5]) = fwhmConfY([3 5])-fwhmConfX(7);
% 	fwhmConfX
% 	fwhmConfY


%	% NOTE: these are the ones used for my ISBI2013 paper! the problem was that this are different coverslips than used with the measured data!!!!!!!
% 	% NOTE: these are calibration results from analyzeRIM (ISBI1) taken at zero depth with combined clens fitting under normal conditions:
% 	% 	-removeFirstSliceInStack = 1
% 	% 	-doCombinedCLensFit = 1
% 	% 	-intensityThreshold = 0.2 % 0.12
% 	% 	-alignmentAlgorithm = 1
% 	% these values are used to simulate mdaSCRIPT2() results if always only the zero depth calibration is taken as the calibration curve. thus, we show that the
% 	% error becomes big in this case
% 	fprintf('\n\nWARNING: resetting calibration results to general calibration at zero depth!\n\n');
% 	fwhmConfX = [262.4401280019837 299.6894840071584 53.0346174334790 0.1322849202155 50.7447745676622 -53.0330558650582 10.1043257693032];
% 	fwhmConfY = [262.4401280019837 299.6894840071584 63.1389432027822 0.1322849202155 60.8491003369654                 0                0];
% 	% align to real vertex
% 	fwhmConfX([3 5]) = fwhmConfX([3 5])-fwhmConfY(5);
% 	fwhmConfY([3 5]) = fwhmConfY([3 5])-fwhmConfY(5);
% 	fwhmConfX
% 	fwhmConfY


% 	% NOTE: these are calibration results from analyze3DScan(!!!!) (ISBI4, so they fit to the measurements) taken at zero depth with combined clens fitting under normal conditions:
% 	% 	-removeFirstSliceInStack = 1
% 	% 	-doCombinedCLensFit = 1
% 	% 	-intensityThreshold = 0.12
% 	% 	-alignmentAlgorithm = 0
% 	% 	-alignmentAlgorithmDg = 1
% 	% these values are used to simulate mdaSCRIPT2() results if always only the zero depth calibration is taken as the calibration curve. thus, we show that the
% 	% error becomes big in this case
% 	fprintf('\n\nWARNING: resetting calibration results to general calibration at zero depth!\n\n');
% 	fwhmConfX = [253.6376891168988   271.9854477540136   25.1719512812952  -0.0106291040719   25.3269658173715  -64.0492221545134   10.7856672239672];
% 	fwhmConfY = [253.6376891168988   271.9854477540136   35.9576185052623  -0.0106291040719   36.1126330413387                   0                   0];
% % 	% alternative configuration with intensityThreshold 0.25 expecting that the fit is better than because we anayway analyse our algorithm only in a small axial environment
% % 	fwhmConfX = [261.2057843573178   303.6514955127300   26.8768787486757   0.1058180835169   24.9949247291980  -63.4553749066945   11.0023118049407];
% % 	fwhmConfY = [261.2057843573178   303.6514955127300   37.8791905536163   0.1058180835169   35.9972365341386                   0                   0];
% 	% align to real vertex
% 	fwhmConfX([3 5]) = fwhmConfX([3 5])-fwhmConfY(5);
% 	fwhmConfY([3 5]) = fwhmConfY([3 5])-fwhmConfY(5);
% 	fwhmConfX
% 	fwhmConfY
	
	
% 	% NOTE: these are calibration results from analyze3DScan(!!!!) (ISBI4, so they fit to the measurements) taken at zero depth with seperate clens fitting under normal conditions:
% 	% 	-removeFirstSliceInStack = 1
% 	% 	-doCombinedCLensFit = 0
% 	% 	-intensityThreshold = 0.12
% 	% 	-alignmentAlgorithm = 0
% 	% 	-alignmentAlgorithmDg = 1
% 	% these values are used to simulate mdaSCRIPT2() results if always only the zero depth calibration is taken as the calibration curve. thus, we show that the
% 	% error becomes big in this case
% 	fprintf('\n\nWARNING: resetting calibration results to general calibration at zero depth!\n\n');
% 	fwhmConfX = [330.9365986973152   357.0408950533418   29.9799434428747   0.1864212598001   26.3145284006302  -68.9794393399248   9.8782970889941];
% 	fwhmConfY = [258.9956099286020   274.9858405247785   33.8109380513862  -0.1607688917119   36.1928254896243                   0                   0];
% % 	% alternative configuration with intensityThreshold 0.25 expecting that the fit is better than because we anayway analyse our algorithm only in a small axial environment
% % 	fwhmConfX = [420.0742775852464   422.7182066263270   39.9986414604667   0.6340573033979   22.4846334431128  -64.4266152238208   13.4686711572730];
% % 	fwhmConfY = [261.5157273280561   303.7756083074182   33.8436777137653  -0.1184345637170   35.9533046003857                   0                   0];
% 	% align to real vertex
% 	fwhmConfX([3 5]) = fwhmConfX([3 5])-fwhmConfY(5);
% 	fwhmConfY([3 5]) = fwhmConfY([3 5])-fwhmConfY(5);
% 	fwhmConfX
% 	fwhmConfY

% 	% NOTE: these are calibration results from analyzeRIM(!!!!) (ISBI4, so they fit to the measurements) taken at zero depth with combined clens fitting under normal conditions:
% 	% 	-removeFirstSliceInStack = 1
% 	% 	-doCombinedCLensFit = 1
% 	% 	-intensityThreshold = 0.12
% 	% 	-alignmentAlgorithm = 1
% 	% these values are used to simulate mdaSCRIPT2() results if always only the zero depth calibration is taken as the calibration curve. thus, we show that the
% 	% error becomes big in this case
% 	fprintf('\n\nWARNING: resetting calibration results to general calibration at zero depth!\n\n');
% 	fwhmConfX = [254.8978764657605   274.7130172848530   25.0060278383613  -0.0134853035829   25.2056776315465  -62.9153903005184   10.8452446554745];
% 	fwhmConfY = [254.8978764657605   274.7130172848530   35.8512724938358  -0.0134853035829   36.0509222870210                   0                   0];
% 	% align to real vertex
% 	fwhmConfX([3 5]) = fwhmConfX([3 5])-fwhmConfY(5);
% 	fwhmConfY([3 5]) = fwhmConfY([3 5])-fwhmConfY(5);
% 	fwhmConfX
% 	fwhmConfY


% 	% NOTE: these are calibration results from analyzeRIM(!!!!) (ISBI4, so they fit to the measurements) taken at zero depth with seperate clens fitting under normal conditions:
% 	% 	-removeFirstSliceInStack = 1
% 	% 	-doCombinedCLensFit = 0
% 	% 	-intensityThreshold = 0.12
% 	% 	-alignmentAlgorithm = 1
% 	% these values are used to simulate mdaSCRIPT2() results if always only the zero depth calibration is taken as the calibration curve. thus, we show that the
% 	% error becomes big in this case
% 	fprintf('\n\nWARNING: resetting calibration results to general calibration at zero depth!\n\n');
% 	fwhmConfX = [331.2743375577679   363.7717395901626   29.5727814619299   0.1830116796340   25.8413898181485  -69.4003525950679   10.1819631258368];
% 	fwhmConfY = [258.2987654707981   273.2932351771887   33.8579744853145  -0.1479254354135   36.0233529439854                   0                   0];
% 	% align to real vertex
% 	fwhmConfX([3 5]) = fwhmConfX([3 5])-fwhmConfY(5);
% 	fwhmConfY([3 5]) = fwhmConfY([3 5])-fwhmConfY(5);
% 	fwhmConfX
% 	fwhmConfY

	

% 	% NOTE: these are calibration results from analyzeRIM (ISBI1) taken at zero depth with separated clens fitting under normal conditions:
% 	% 	-removeFirstSliceInStack = 1
% 	% 	-doCombinedCLensFit = 0
% 	% 	-intensityThreshold = 0.2 % 0.12
% 	% 	-alignmentAlgorithm = 1
% 	% these values are used to simulate mdaSCRIPT2() results if always only the zero depth calibration is taken as the calibration curve. thus, we show that the
% 	% error becomes big in this case
% 	fprintf('\n\nWARNING: resetting calibration results to general calibration at zero depth!\n\n');
% 	fwhmConfX = [470.4750436931768   381.2323597493635   69.9001518949993   0.9014136470012   49.5145863602762  -62.6312566705284   11.3390643179306];
% 	fwhmConfY = [259.4039980849602   298.9263959788606   59.7490375687020  -0.0639595634401   60.8536506782068                   0                   0];
% 	% remove offsets because we fitted separately
% 	fwhmConfX([6 7]) = 0;
% 	% align to real vertex
% 	fwhmConfX([3 5]) = fwhmConfX([3 5])-fwhmConfY(5);
% 	fwhmConfY([3 5]) = fwhmConfY([3 5])-fwhmConfY(5);
% 	fwhmConfX
% 	fwhmConfY

% 	% NOTE: these are calibration results for WIDE-FIELD MODE from analyzeRIM taken at zero depth with combined clens fitting under normal conditions:
% 	% 	-removeFirstSliceInStack = 1
% 	% 	-doCombinedCLensFit = 1
% 	% 	-intensityThreshold = 0.7 % 0.12
% 	% 	-alignmentAlgorithm = 1
% 	% these values are used to simulate mdaSCRIPT2() results if always only the zero depth calibration is taken as the calibration curve. thus, we show that the
% 	% error becomes big in this case
% 	fprintf('\n\nWARNING: resetting calibration results to general calibration at zero depth!\n\n');
% 	fwhmConfX = [295.1786889974718   286.8598235868252   23.5054353733199  -0.1205860289124   25.1979178753381  -11.5451770551526   11.0971654990259];
% 	fwhmConfY = [295.1786889974718   286.8598235868252   34.6026008723458  -0.1205860289124   36.2950833743640                   0                   0];
% 	% align to real vertex
% 	fwhmConfX([3 5]) = fwhmConfX([3 5])-fwhmConfY(5);
% 	fwhmConfY([3 5]) = fwhmConfY([3 5])-fwhmConfY(5);
% 	fwhmConfX
% 	fwhmConfY

% 	% NOTE: these are calibration results from analyzeRIM (ISBI1) taken at zero depth with combined clens fitting AND forced symmetrical parabola under normal conditions:
% 	% 	-removeFirstSliceInStack = 1
% 	% 	-doCombinedCLensFit = 1
% 	% 	-intensityThreshold = 0.2 % 0.12
% 	% 	-alignmentAlgorithm = 1
% 	% these values are used to simulate mdaSCRIPT2() results if always only the zero depth calibration is taken as the calibration curve. thus, we show that the
% 	% error becomes big in this case
% 	fprintf('\n\nWARNING: resetting calibration results to general calibration at zero depth!\n\n');
% 	fwhmConfX = [255.9039861884350   288.9567134606074   50.0333298168922   0.0000000099797   50.0333296540846  -55.5042125345272   10.1166572668664];
% 	fwhmConfY = [255.9039861884350   288.9567134606074   60.1499870837586   0.0000000099797   60.1499869209510                   0                   0];
% 	% align to real vertex
% 	fwhmConfX([3 5]) = fwhmConfX([3 5])-fwhmConfY(5);
% 	fwhmConfY([3 5]) = fwhmConfY([3 5])-fwhmConfY(5);
% 	fwhmConfX
% 	fwhmConfY

end

%% compute the full 3D information for each particle

if computeZPosition
	
% compute lookup table for the synthetic calibration function

% resolution of the lookup table to compute the inverse of the difference function for the cylindrical Lens
% here i want 1 nm resolution
axialLTResInNm = 1;
axialInverseResolution = axialLTResInNm/dv(3);
axialSize = 2000;%100;

% compute the difference function of the biased fwhms with high resolution, which is later used
% as a lookup table to compute the inverse function
% NOTE: that the difference function is only UNAMBIGUOUS if the two fwhm curves cross only ONCE and this only happens if they have the same configuration and
% there are ONLY shifts!!!!!!!!!!!!!!!!!!!!!!!!!!!
axialGrid = [-axialSize:axialInverseResolution:axialSize];
axialLookupTable = biasedFwhmDifferenceFcn1DDiscrete(fwhmConfX, fwhmConfY, axialGrid, dv(3), doCombinedCLensFit);

% % plot reference and estimated calibration curves and difference functions
% figure;
% hold on;
% title('curves');
% % if doCombinedCLensFit
% % 	% reference curves
% % 	plot(axialGrid, biasedFwhm1DDiscrete(fwhmRefConfX, axialGrid, dv(3))-fwhmRefConfX(6), 'r');
% % 	plot(axialGrid, biasedFwhm1DDiscrete(fwhmRefConfY, axialGrid, dv(3)), 'r');
% % % 	plot(axialGrid, biasedFwhmDifferenceFcn1DDiscrete(fwhmRefConfX, fwhmRefConfY, axialGrid, dv(3)), 'b');
% % 	% reference coordinates
% % 	plot(fwhmRefCurveCoord(1:3, 1), fwhmRefCurveCoord(1:3, 2), 'rs', 'LineStyle', 'none');
% % else
% % 	% reference curves
% % 	plot(axialGrid, biasedFwhm1DDiscrete(fwhmRefConfX, axialGrid, dv(3))+fwhmRefConfX(6), 'r');
% % 	plot(axialGrid, biasedFwhm1DDiscrete(fwhmRefConfY, axialGrid, dv(3)), 'r');
% % % 	plot(axialGrid, biasedFwhmDifferenceFcn1DDiscrete(fwhmRefConfX, fwhmRefConfY, axialGrid, dv(3)), 'b');
% % 	% reference coordinates
% % 	plot(fwhmRefCurveCoord(1:3, 1), fwhmRefCurveCoord(1:3, 2), 'rs', 'LineStyle', 'none');
% % end
% % estimated curves
% % plot(axialGrid, biasedFwhm1DDiscrete(fwhmConfX, axialGrid, dv(3))-fwhmConfX(6), 'k');
% plot(axialGrid, biasedFwhm1DDiscrete(fwhmConfX, axialGrid, dv(3)), 'k');
% plot(axialGrid, biasedFwhm1DDiscrete(fwhmConfY, axialGrid, dv(3)), 'k');
% % plot(axialGrid, axialLookupTable, 'g');
% % estimated coordinates
% plot(fwhmCurveCoord(1:3, 1), fwhmCurveCoord(1:3, 2), 'ks', 'LineStyle', 'none');
% hold off;

% loop over frames
for fNum = 1:size(p, 1)%numel(p)%size(im, 3)
	
	% loop over particles
	for pNum = 1:size(p{fNum}, 1)

		% if there is no particle
		if pNum == 0
			% just continue
			continue;
		end
		
		%% compute z values using a lookup table computed from the calibration curves
		% NOTE: that the difference function is only UNAMBIGUOUS if the two fwhm curves cross only ONCE and this only happens if they have the same configuration and
		% there are ONLY shifts!!!!!!!!!!!!!!!!!!!!!!!!!!!
		p{fNum}(pNum, 10) = inverseBiasedFwhmDifferenceFcn1D(p{fNum}(pNum, 3)*dv(1)-p{fNum}(pNum, 4)*dv(2), axialGrid, axialLookupTable).*dv(3);
% 		p{fNum}(pNum, 10) = inverseBiasedFwhmDifferenceFcn1D(((p{fNum}(pNum, 3)*dv(1))+fwhmConfX(6))-p{fNum}(pNum, 4)*dv(2), axialGrid, axialLookupTable).*dv(3);
		
		%% compute z values by fusion of the xy fwhm values

		% compute z values from xy fwhm
		% NOTE: because the parameters of the fwhm in x direction do not include the fwhm offset we simply add that offset to the obtained y coordinate, which
		% has the same effect
		if doCombinedCLensFit
			zx = inverseBiasedFwhm1D(fwhmConfX, (p{fNum}(pNum, 3)*dv(1))+fwhmConfX(6), dv(3));
		else
			zx = inverseBiasedFwhm1D(fwhmConfX, (p{fNum}(pNum, 3)*dv(1)), dv(3));
		end
% 		zx = inverseBiasedFwhm1D(fwhmConfX, p{fNum}(pNum, 3)*dv(1)+fwhmConfX(6), dv(3));
		zy = inverseBiasedFwhm1D(fwhmConfY, p{fNum}(pNum, 4)*dv(2), dv(3));

		% merge the two z values that are closest
		% NOTE: this is acceptable as the two fwhm curves are well separated. however it creates uncertainty in the vicinity of the vertex points of the fwhm
		% because here the two values from one curve are too close to be clearly connected to the point at the other fwhm curve
		% NOTE: as we have only 4 values i hack it directly
		[null closestValuesIdx] = min(abs([zx(1)-zy(1) zx(1)-zy(2) zx(2)-zy(1) zx(2)-zy(2)]));
		
		switch closestValuesIdx
			case 1
				p{fNum}(pNum, 11) = mean([zx(1) zy(1)]);
			case 2
				p{fNum}(pNum, 11) = mean([zx(1) zy(2)]);
			case 3
				p{fNum}(pNum, 11) = mean([zx(2) zy(1)]);		
			case 4
				p{fNum}(pNum, 11) = mean([zx(2) zy(2)]);
		end
		
		%% compute z value by finding the minimal distance between the calibration curves and the measured fwhm values in the square root space (see suppl. huang2008)
		% NOTE: the search for z is initialized with the result from the z value created via fusion
		p{fNum}(pNum, 11) = findAxialPositionInSquareRootSpaceOfMeasuredFwhms(p{fNum}(pNum, 11), p{fNum}(pNum, 3)*dv(1), p{fNum}(pNum, 4)*dv(2), fwhmConfX, fwhmConfY, dv, axialRange, 0.1, doCombinedCLensFit);

% 		% choose direction relative to the particle by analysis of the difference between xy fwhm due to the lens
% 		% NOTE: this is NONSENSE as the crossing of the curves has nothing to do with the center of the particle. indeed the center is at the non-shifted fwhm
% 		% curve
% 		if p{fNum}(pNum, 3)*dv(1) - p{fNum}(pNum, 4)*dv(2) >= 0
% 			p{fNum}(pNum, 11) = mean([zx(2) zy(2)]);
% 		else
% 			p{fNum}(pNum, 11) = mean([zx(1) zy(1)]);
% 		end

	end
end

end

%%% connect to traces using the 3D information, data fusion, blinking, motion model and fit quality

end

function y = equationSystem(par)

global fwhmCurveCoord;
sign = -1;

bevor man das hier benutzt muss man nochmal auf die axialResolution achten die ich hier manchmal mit 20 drin hab

% y = [par(1)*sqrt(1+(fwhmCurveCoord(1,1)/par(2))^2)+par(3)*fwhmCurveCoord(1,1) - fwhmCurveCoord(1,2);...
% 	par(1)*sqrt(1+(fwhmCurveCoord(2,1)/par(2))^2)+par(3)*fwhmCurveCoord(2,1) - fwhmCurveCoord(2,2);...
% 	par(1)*sqrt(1+(fwhmCurveCoord(3,1)/par(2))^2)+par(3)*fwhmCurveCoord(3,1) - fwhmCurveCoord(3,2)];

y = [(par(1)/(par(2)^2*sqrt(1+(fwhmCurveCoord(1,1)/par(2))^2)))*fwhmCurveCoord(1,1)+par(3);...
	par(1)*sqrt(1+(fwhmCurveCoord(2,1)/par(2))^2)+par(3)*fwhmCurveCoord(2,1) - fwhmCurveCoord(2,2);...
	par(1)*sqrt(1+(fwhmCurveCoord(3,1)/par(2))^2)+par(3)*fwhmCurveCoord(3,1) - fwhmCurveCoord(3,2)];

% y = [par(1)*sqrt(1+((fwhmCurveCoord(1,1)+sign*sqrt((par(3)^2*par(2)^4)/(par(1)^2-par(3)^2*par(2)^2))/20)/par(2))^2)+par(3)*(fwhmCurveCoord(1,1)+sign*sqrt((par(3)^2*par(2)^4)/(par(1)^2-par(3)^2*par(2)^2))/20) - fwhmCurveCoord(1,2);...
% 	par(1)*sqrt(1+((fwhmCurveCoord(2,1)+sign*sqrt((par(3)^2*par(2)^4)/(par(1)^2-par(3)^2*par(2)^2))/20)/par(2))^2)+par(3)*(fwhmCurveCoord(2,1)+sign*sqrt((par(3)^2*par(2)^4)/(par(1)^2-par(3)^2*par(2)^2))/20) - fwhmCurveCoord(2,2);...
% 	par(1)*sqrt(1+((fwhmCurveCoord(3,1)+sign*sqrt((par(3)^2*par(2)^4)/(par(1)^2-par(3)^2*par(2)^2))/20)/par(2))^2)+par(3)*(fwhmCurveCoord(3,1)+sign*sqrt((par(3)^2*par(2)^4)/(par(1)^2-par(3)^2*par(2)^2))/20) - fwhmCurveCoord(3,2)];

% y = [(par(1)/(par(2)^2*sqrt(1+((fwhmCurveCoord(1,1)+sign*sqrt((par(3)^2*par(2)^4)/(par(1)^2-par(3)^2*par(2)^2))/20)/par(2))^2)))*(fwhmCurveCoord(1,1)+sign*sqrt((par(3)^2*par(2)^4)/(par(1)^2-par(3)^2*par(2)^2))/20)+par(3);...
% 	par(1)*sqrt(1+((fwhmCurveCoord(2,1)+sign*sqrt((par(3)^2*par(2)^4)/(par(1)^2-par(3)^2*par(2)^2))/20)/par(2))^2)+par(3)*(fwhmCurveCoord(2,1)+sign*sqrt((par(3)^2*par(2)^4)/(par(1)^2-par(3)^2*par(2)^2))/20) - fwhmCurveCoord(2,2);...
% 	par(1)*sqrt(1+((fwhmCurveCoord(3,1)+sign*sqrt((par(3)^2*par(2)^4)/(par(1)^2-par(3)^2*par(2)^2))/20)/par(2))^2)+par(3)*(fwhmCurveCoord(3,1)+sign*sqrt((par(3)^2*par(2)^4)/(par(1)^2-par(3)^2*par(2)^2))/20) - fwhmCurveCoord(3,2)];

end

%% function y = biasedFwhm1DFunc (init, x)
%
% frontend to biasedFwhm1DDiscrete() to be able use it with lsqcurvefit()
% NOTE: this function assumes an integer grid thus all inputs have to be defined as integer as well
% NOTE: after schuetz2000
%
% parameters:
% -init		- vector of parameters defining the initial guess of the parabola:
%  		   -(1) - fwhm in focus (thus the resolution AND the angular point of the parabol)
%          -(2) - depth of focus
%		   -(3) - shifting of the parabol at the z-scale
%		   -(4) - gradient (m) of the straight line (mx+n) that biases the fwhm (no offset n)
%		   -..	- further parameters are up to you but are ignored here
% -x      - x values (axial z)
%
% returns:
% -y		- y values of the parabola (fwhm off focus)

function y = biasedFwhm1DFunc (init, x)

global axialRes;
% global vertexPointCoord;

% y = biasedFwhm1DDiscrete(init, x, axialRes);
% y = biasedFwhm1DDiscrete([init(1) init(2) vertexPointCoord init(3)], x, axialRes);
y = biasedFwhm1DDiscrete([init(1) init(2) sqrt((init(3)^2*init(2)^4)/(init(1)^2-init(3)^2*init(2)^2))/axialRes init(3)], x, axialRes);

end