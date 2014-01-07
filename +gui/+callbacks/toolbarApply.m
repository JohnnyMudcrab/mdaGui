function [] = toolbarApply(this)

  %% set status to busy
  this.changeStatus('statusMain', 'Busy...');
  
  %% get tree
  hTree = this.getHandle('treeMain');
  
  selectedNodes = hTree.getSelectedNodes();
  
  % non selected
  if isempty(selectedNodes)
    msgbox('No Image selected','No Image','Help')
    return
  end
  
  % more then one selected
  if numel(selectedNodes) > 1 || ~selectedNodes{1}.isLeaf
    msgbox('Too many Images selected','Too many Images','Help')
    return
  end
  
  %% save current Node
  hMda = getappdata(0,'hMda');
  hMda.currentNode = selectedNodes{1};
  
  
  %% get file location
  file = selectedNodes{1}.handle.UserData.string;
  path = selectedNodes{1}.getParent.handle.UserData.string;
  
  
  %% show picture in axes
  info = imfinfo([path file]); 
  img = imread([path file],'Index',1,'Info',info); 
  
  axes(this.getHandle('axesTrack'))
  imshow(double(img) / 255) 
  
  axes(this.getHandle('axesLocate'))
  imshow(double(img) / 255) 
  
  axes(this.getHandle('axesInfo'))
  imshow(double(img) / 255)
    
  
  %% update info panel
  set(this.getHandle('textInfoName2'), 'String', [path file])
  set(this.getHandle('textInfoResolution2'), 'String', [num2str(info(1, 1).Width) ' x ' num2str(info(1, 1).Height)])
  set(this.getHandle('textInfoNumber2'), 'String', num2str(numel(info)))
  set(this.getHandle('textInfoBit2'), 'String', [num2str(info(1,1).BitDepth) ' Bit'])
  
  
  %% update state
  if ~isempty(selectedNodes{1}.handle.UserData.locate)
    set(this.getHandle('textInfoStateLocate'), 'String', 'Located')
  else
    set(this.getHandle('textInfoStateLocate'), 'String', 'Not Located Yet')
  end
  
  if ~isempty(selectedNodes{1}.handle.UserData.locate)
    set(this.getHandle('textInfoStateTrack'), 'String', 'Tracked')
  else
    set(this.getHandle('textInfoStateTrack'), 'String', 'Not Tracked Yet')
  end
  
  if ~isempty(selectedNodes{1}.handle.UserData.locate)
    set(this.getHandle('textInfoStateCalibrate'), 'String', 'Calibrated')
  else
    set(this.getHandle('textInfoStateCalibrate'), 'String', 'Not Calibrated Yet')
  end
  
  
  %% enable slider info
  set(this.getHandle('textInfoFrame'), 'String', ['1 / ' num2str(numel(info))])
  set(this.getHandle('sliderInfo'), 'Enable', 'On');
  sliderStep = [1 10] / (numel(info) - 1);
  set(this.getHandle('sliderInfo'), 'Min', 1, 'Max', numel(info), 'SliderStep', sliderStep, 'Value', 1)
  
  
  %% enable locate
  set(this.getHandle('textLocateFrame'), 'String', ['1 / ' num2str(numel(info))])
  set(this.getHandle('sliderLocate'), 'Enable', 'On');
  
  sliderStep = [1 10] / (numel(info) - 1);
  set(this.getHandle('sliderLocate'), 'Min', 1, 'Max', numel(info), 'SliderStep', sliderStep, 'Value', 1)
  set(this.getHandle('buttonLocateApply'), 'Enable', 'on')
  set(this.getHandle('buttonLocateApplyAll'), 'Enable', 'on')
  
  %% enable track
  set(this.getHandle('textTrackFrame'), 'String', ['1 / ' num2str(numel(info))])
  set(this.getHandle('sliderTrack'), 'Enable', 'On');
  
  sliderStep = [1 10] / (numel(info) - 1);
  set(this.getHandle('sliderTrack'), 'Min', 1, 'Max', numel(info), 'SliderStep', sliderStep, 'Value', 1)
  set(this.getHandle('buttonTrack'), 'Enable', 'on')
    
  
  %% set status to ready
  this.changeStatus('statusMain', 'Ready...');

end

