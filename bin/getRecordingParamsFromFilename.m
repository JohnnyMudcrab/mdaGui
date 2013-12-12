function [dv EMGain frameRate binning lengthTS depthInMedia zRange stimDelay zScan isStepScan preMagnification] = getRecordingParamsFromFilename(filename, setup, debug)

dv = [];
EMGain = [];
frameRate = [];
binning = [];
lengthTS = [];
depthInMedia = [];
zRange = [];
stimDelay = [];
zScan= [];

%% extract dv parameter
filename2 = filename(strfind(filename, 'OIL')+3:end);
if length(filename2) == 0
	filename2 = filename(strfind(filename, 'H2O')+3:end);

	if length(filename2) == 0
		filename2 = filename(strfind(filename, 'TIRF')+4:end);
	end
end
filename3 = filename2(1:strfind(filename2, '_')-1);
magnification = str2num(filename3);

% extract pre magnification lens
filename2 = filename(strfind(filename, 'MagLens')+7:end);
if length(filename2) > 0
	filename3 = filename2(1:strfind(filename2, '_')-1);
	preMagnification = str2num(filename3);
else
	preMagnification = 1;
end

% compute according dv
% if it is my confocal setup:
if setup == 1
	dv(1) = round(1000/((magnification*preMagnification*1.2)/16));
	dv(2) = dv(1);
% martins setup
elseif setup == 2
	dv(1) = round(1000/((magnification*preMagnification)/16));
	dv(2) = dv(1);
% if it is my confocal setup BUT with andors NEO camera:
elseif setup == 3
	dv(1) = round(1000/((magnification*preMagnification*1.2)/6.5));
	dv(2) = dv(1);
% Tanimas widefield with small pixel camera
elseif setup == 4
	dv(1) = round(1000/((magnification*preMagnification)/8));
	dv(2) = dv(1);
else
	fprintf('ERROR: getRecordingParamsFromFilename: unknown setup indicator\n')
end
% if we have a 3D file
filename2 = filename(strfind(filename, 'ZSCAN')+5:end);
if length(filename2) > 0
	filename3 = filename2(1:strfind(filename2, '_')-1);
	dv(3) = str2num(filename3);
end

%% extract EMGain parameter
filename2 = filename(strfind(filename, 'EM')+2:end);
filename3 = filename2(1:strfind(filename2, '_')-1);
EMGain = str2num(filename3);

%% extract frame rate parameter
% filename2 = filename(strfind(filename, 'FR')+2:end);
% filename3 = filename2(1:strfind(filename2, '_')-1);
% frameRate = str2num(filename3);
filename2 = filename(strfind(filename, 'ET')+2:end);
filename3 = filename2(1:strfind(filename2, '_')-1);
frameRate = 1000/str2num(filename3);

%% extract binning parameter
filename2 = filename(strfind(filename, 'BIN')+3:end);
filename3 = filename2(1:strfind(filename2, '_')-1);
binning = str2num(filename3);

%% extract length parameter for TS images
filename2 = filename(strfind(filename, 'TS')+2:end);
if length(filename2) > 0
	filename3 = filename2(1:strfind(filename2, '_')-1);
	lengthTS = str2num(filename3);
else
	lengthTS = [];
end

%% extract depth in media parameter
filename2 = filename(strfind(filename, 'DIST')+4:end);
filename3 = filename2(1:strfind(filename2, '_')-1);
depthInMedia = str2num(filename3);

%% extract zRange parameter
filename2 = filename(strfind(filename, 'ZRANGE')+6:end);
filename3 = filename2(1:strfind(filename2, '_')-1);
zRange = str2num(filename3);

%% extract stimulation delay parameter
filename2 = filename(strfind(filename, 'STIM')+4:end);
filename3 = filename2(1:strfind(filename2, '_')-1);
stimDelay = str2num(filename3);

%% extract zScan parameter
filename2 = filename(strfind(filename, 'ZSCAN')+5:end);
filename3 = filename2(1:strfind(filename2, '_')-1);
zScan = str2num(filename3);

%% extract isStepScan parameter
filename2 = filename(strfind(filename, 'StepScan')+8:end);
if isempty(filename2)
	isStepScan = 0;
else
	isStepScan = 1;
end

if debug
	if length(dv)==2
		fprintf('took file with dv: [%d %d] EMGain: %d frame rate: %d binning: %d lengthTS: %d depthInMedia: %d zRange: %d stimDelay: %d zScan: %d isStepScan: %d preMagnification: %d\n', dv(1), dv(2), EMGain, frameRate, binning, lengthTS, depthInMedia, zRange, stimDelay, zScan, isStepScan, preMagnification);
	elseif length(dv)==3
		fprintf('took file with dv: [%d %d %d] EMGain: %d frame rate: %d binning: %d depthInMedia: %d zRange: %d stimDelay: %d zScan: %d isStepScan: %d preMagnification: %d\n', dv(1), dv(2), dv(3), EMGain, frameRate, binning, depthInMedia, zRange, stimDelay, zScan, isStepScan, preMagnification);
	end
end

end