% im - image
% iDiff - intensity difference required between top and brim region
% dTop, dBrim - diameter  (integer number)
function imTH = topHatFilter(im, iDiff, dTop, dBrim, createBinaryIm, shrinkBinaryIm2Point)

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
end

% loop over image
for x = ceil(dBrim/2):size(im, 2)-thHalf
	for y = ceil(dBrim/2):size(im, 1)-thHalf
		
		% cut filter region from image
		imCut = im(y-thHalf:y+thHalf, x-thHalf:x+thHalf);

		if max(imCut(thTop))-max(imCut(thBrim)) >= iDiff 
			imTH(y, x) = im(y, x);
		end
	end
end

if debug
	figure,imshow(imTH, []);
end

% binearize image if requested
if createBinaryIm

	imTH(imTH>0) = 1;

	if debug
		figure,imshow(imTH, []);
	end

	% shrink to point if requested
	if shrinkBinaryIm2Point
		
		imTH = bwmorph(imTH, 'Shrink', Inf);
		
		if debug
			figure,imshow(imTH, []);
		end

	end
end

end