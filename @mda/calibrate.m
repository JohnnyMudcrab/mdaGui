function calibrate(this)
  
  hTree = this.gui.getHandle('treeMain');


  [selectedNodes, ~, ~] = hTree.getSelectedNodes();


% minTrack = str2double(get(handles.edit_minimal_track,'String')); % get minTarck from gui
%  
% voxelSize = [str2double(get(handles.edit_voxelsize_x,'String')),...
%              str2double(get(handles.edit_voxelsize_y,'String')),1];
%            
% fwhm_offset = str2double(get(handles.edit_fwhm_offset,'String'));
% lens_shift = str2double(get(handles.edit_lens_shift,'String'));

%% Make List of images
[images,nodes, nodeIndex] = mda_makeList(selectedNodes);


%% loop over images

imageCount = size(images,1);

for nodeIterator = 1:1:imageCount
    nodeIterator
  currentNode = nodes{nodeIndex(nodeIterator,1),1};
  
  data = currentNode.handle.UserData;


  %% create p_container and delete short path

  p_container = data.track.trajectories;
  trackingParams = create_trackingParam;

  % default
  fwhmConfX = [0    0   0   0    0   -data.conf.fwhm_offset  data.conf.lens_shift];
  fwhmConfY = [0    0   0   0    0       0              0];

  n = size(p_container,1);

%   for i = n:-1:1
%     if p_container{i,3} < data.conf.minTrack
%       p_container(i,:) = [];
%     end
%   end

  %% loop over p

  combCalibration = get(handles.checkbox_combCalibration,'Value');


  if (combCalibration)

    % combine tracks
    temp  = p_container(:,1);
    combined_p = {};
    for i = 1:size(temp,1)
      if i==9
        disp('test')
      end
      combined_p = [combined_p; temp{i,1}(:,1:2); {[] []}];
    end

    clear temp
	

% 	% combine tracks but only include tracks that show reasonable lateral movement
% 	fprintf('\nNOTE: combining tracks but only these having reasonable lateral movement! \n');
% 	combined_p = {};
% 
%     temp  = p_container(:,1);
% 	% loop over tracks
%     for i = 1:size(temp, 1)
% 		% compute the lateral movement of this track
% 		minx = 100000000;
% 		maxx = 0;
% 		miny = 100000000;
% 		maxy = 0;
% 
% 		% get track
% 		null = temp{i, 1};
% 		minMaxHasBeenAdjustedFlag = 0;
% 		
% 		% find min and max values of the lateral movement of this track
% 		for n  = 1:size(null, 1)
% 			if ~isempty(null{n})
% 				minMaxHasBeenAdjustedFlag = 1;
% 				if minx>null{n}(7);	minx = null{n}(7); end;
% 				if miny>null{n}(8);	miny = null{n}(8); end;
% 				if maxx<null{n}(7);	maxx = null{n}(7); end;
% 				if maxy<null{n}(8);	maxy = null{n}(8); end;
% 			end
% 		end
% 		
% 		if minMaxHasBeenAdjustedFlag
% 			nullLength = sqrt((maxx-minx)^2 + (maxy-miny)^2);
% 		else
% 			nullLength = 0;
% 		end
% 		
% 		% if it travelled long enough to be considered as being not stationary
% 		if nullLength>4
% 			% add it for calibration
% 			combined_p = [combined_p; temp{i,1}; {[] []}];
% 		end
%     end
% 
%     clear temp

% 	% combine tracks but only include tracks that show reasonable lateral movement
% 	% NOTE: this is the old version for NON INDEXED tracks (NoCaliHack)
% 	fprintf('\nNOTE: combining tracks but only these having reasonable lateral movement (NoCaliHack)!\n');
% 	combined_p = {};
% 
% 	% loop over tracks
% 	for i = 1:size(p_container,1)
% 		minx = 100000000;
% 		maxx = 0;
% 		miny = 100000000;
% 		maxy = 0;
% 
% 		null = p_container{i, 1};
% 		minMaxHasBeenAdjustedFlag = 0;
% 		
% 		% find min and max values of the lateral movement of this track
% 		for n  = 1:size(null, 1)
% 			if ~isempty(null{n})
% 				minMaxHasBeenAdjustedFlag = 1;
% 				if minx>null{n}(7);	minx = null{n}(7); end;
% 				if miny>null{n}(8);	miny = null{n}(8); end;
% 				if maxx<null{n}(7);	maxx = null{n}(7); end;
% 				if maxy<null{n}(8);	maxy = null{n}(8); end;
% 			end
% 		end
% 		
% 		if minMaxHasBeenAdjustedFlag
% 			nullLength = sqrt((maxx-minx)^2 + (maxy-miny)^2);
% 		else
% 			nullLength = 0;
% 		end
% 		
% 		% if it travelled long enough to be considered as being not stationary
% 		if nullLength>3
% 			% add it for calibration
% 			combined_p = [combined_p; null];
% 		end
% 	end
% 	clear null;

	if isempty(combined_p)
		fwhmConfX = [];
		fwhmConfY = [];
		success = -1;
		control = {};
	else
		% do calibration
		[~, fwhmConfX, fwhmConfY, ~, ~, success, control] = mda2([], data.conf.voxelSize, trackingParams, 0, 64, combined_p,fwhmConfX, fwhmConfY,[-500 - data.conf.lens_shift, 500] ./ data.conf.voxelSize(3) ,1, 1, 0, 1,0);
	end
  end

    % do calibration 

  n = size(p_container,1);
  p = cell(n,1);



  for i = 1:1:n

    if(combCalibration && success >= 0)
      [p_current, fwhmConfX, fwhmConfY, ~, ~, ~, ~] = mda2([], data.conf.voxelSize, trackingParams, 0, 64, p_container{i},fwhmConfX, fwhmConfY,[-500 - data.conf.lens_shift, 500] ./ data.conf.voxelSize(3) ,1, 0 ,0 ,1,1);
    else
      if (combCalibration && success == -1)
        p_current = p_container{i};
      end
    end

    if(~combCalibration)
      [p_current, fwhmConfX, fwhmConfY, ~, ~, ~, ~] = mda2([], data.conf.voxelSize, trackingParams, 0, 64, p_container{i},fwhmConfX, fwhmConfY,[-500 - data.conf.lens_shift, 500] ./ data.conf.voxelSize(3) ,1, 1, 0, 1,1);
    end

    p{i,1} = p_current;
    p{i,2} = fwhmConfX;
    p{i,3} = fwhmConfY;
    p{i,4} = p_container{i,2};
    p{i,5} = p_container{i,3};
	if i == 1
		p{i,6} = control;
	else
		p{i,6} = {};
	end
	
    pos = get(waitbarOverall,'position'); 
    pos(3) = i / n * 155; 
    set(waitbarOverall,'position',pos,'string',sprintf('%.0f%%',i /n * 100)) % update waitbar

    drawnow

  end



  data.calibrate = p;
  currentNode.handle.UserData = data;
  
end


%% clear up

delete(waitbarCurrent) % delete waitbar
delete(waitbarOverall)


end

