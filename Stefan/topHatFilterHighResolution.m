% NOTE: this is an advanced version of the normal top hat filter. differences are:
%	-the mean of the top and brim regions is compared instead of the max value of each region, which decreases the impact of high frequent noise (salt'n'pepper)
%	-instead of binearizing and shrinking we take all the local maxima of each region to binearize the image. this has the following advantage:
% 		->even nearby particles can be distinguished because otherwise two regions would become one region during shrinking as they overlap und thus usually one
% 		instead of two peaks would be found. so the limit of resolution is not set by the size of the top anymore.
% im - image
% iDiff - intensity difference required between top and brim region
% dTop, dBrim - diameter  (integer number)
function imTH = topHatFilterHighResolution(im, iDiff, dTop, dBrim, findLocalMaxima)

debug = 0;

if ~mod(dTop, 2) || ~mod(dBrim, 2)
	fprintf('topHatFilter: diameters must be odd!');
	imTH = [];
	return;
end

% compute filter
thBrim = gaussian2D([ceil(dBrim/2), ceil(dBrim/2), dBrim, dBrim, 1, 0], [1 1], dBrim, dBrim);
thBrim(thBrim<=0.5) = 0;
thBrim(thBrim>0) = 1;
thTop = gaussian2D([ceil(dBrim/2), ceil(dBrim/2), dTop, dTop, 1, 0], [1 1], dBrim, dBrim);
thTop(thTop<=0.5) = 0;
thTop(thTop>0) = 1;
thBrim = thBrim-thTop;
thTop = cast(thTop, 'logical');
thBrim = cast(thBrim, 'logical');

% init filtered image
imTH = zeros(size(im));

% compute discretized halfSize of the filter
thHalf = floor(dBrim/2);

if debug
	thBrim
	thTop
end

% loop over image
for x = ceil(dBrim/2):size(im, 2)-thHalf
	for y = ceil(dBrim/2):size(im, 1)-thHalf
		
		% cut filter region from image
		imCut = im(y-thHalf:y+thHalf, x-thHalf:x+thHalf);

% 		if max(imCut(thTop))-max(imCut(thBrim)) >= iDiff 
 		if mean(imCut(thTop))-mean(imCut(thBrim)) >= iDiff 
			imTH(y, x) = im(y, x);
		end
	end
end

if debug
	figure,imshow(imTH, []);
end

% if finding all local maxima in the regions und, thus, binearizing the image is requested
if findLocalMaxima
	
	if sum(imTH(:)) > 0
		imTH = imregionalmax(imTH, 8);
	end
	
	if debug
		figure,imshow(imTH, []);
	end
end

end