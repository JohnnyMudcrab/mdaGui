function apply(this)

  %% set status to busy
  this.gui.changeStatus('statusMain', 'Busy...');
  
  %% get tree
  hTree = this.gui.getHandle('treeMain');
  
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
  this.currentNode = selectedNodes{1};
  
  
  %% get file location
  file = selectedNodes{1}.handle.UserData.string;
  path = selectedNodes{1}.getParent.handle.UserData.string;
  
  
  %% show picture in axes
  info = imfinfo([path file]); 
  img = imread([path file],'Index',1,'Info',info); 
  
  axes(this.gui.getHandle('axesTrack'))
  imshow(double(img) / 255) 
  
  axes(this.gui.getHandle('axesLocate'))
  imshow(double(img) / 255) 
  
  axes(this.gui.getHandle('axesInfo'))
  imshow(double(img) / 255)
    
  
  %% update info panel
  set(this.gui.getHandle('textInfoName2'), 'String', [path file])
  set(this.gui.getHandle('textInfoResolution2'), 'String', [num2str(info(1, 1).Width) ' x ' num2str(info(1, 1).Height)])
  set(this.gui.getHandle('textInfoNumber2'), 'String', num2str(numel(info)))
  set(this.gui.getHandle('textInfoBit2'), 'String', [num2str(info(1,1).BitDepth) ' Bit'])
  
  
  %% update state
  if ~isempty(selectedNodes{1}.handle.UserData.locate)
    set(this.gui.getHandle('textInfoStateLocate'), 'String', 'Located')
  else
    set(this.gui.getHandle('textInfoStateLocate'), 'String', 'Not Located Yet')
  end
  
  if ~isempty(selectedNodes{1}.handle.UserData.locate)
    set(this.gui.getHandle('textInfoStateTrack'), 'String', 'Tracked')
  else
    set(this.gui.getHandle('textInfoStateTrack'), 'String', 'Not Tracked Yet')
  end
  
  if ~isempty(selectedNodes{1}.handle.UserData.locate)
    set(this.gui.getHandle('textInfoStateCalibrate'), 'String', 'Calibrated')
  else
    set(this.gui.getHandle('textInfoStateCalibrate'), 'String', 'Not Calibrated Yet')
  end
  
  
  %% enable slider info
  set(this.gui.getHandle('textInfoFrame'), 'String', ['1 / ' num2str(numel(info))])
  set(this.gui.getHandle('sliderInfo'), 'Enable', 'On');
  sliderStep = [1 10] / (numel(info) - 1);
  set(this.gui.getHandle('sliderInfo'), 'Min', 1, 'Max', numel(info), 'SliderStep', sliderStep, 'Value', 1)
  
  
  %% enable locate
  set(this.gui.getHandle('textLocateFrame'), 'String', ['1 / ' num2str(numel(info))])
  set(this.gui.getHandle('sliderLocate'), 'Enable', 'On');
  
  sliderStep = [1 10] / (numel(info) - 1);
  set(this.gui.getHandle('sliderLocate'), 'Min', 1, 'Max', numel(info), 'SliderStep', sliderStep, 'Value', 1)
  set(this.gui.getHandle('buttonLocateApply'), 'Enable', 'on')
  set(this.gui.getHandle('buttonLocateApplyAll'), 'Enable', 'on')
  
  %% enable track
  set(this.gui.getHandle('textTrackFrame'), 'String', ['1 / ' num2str(numel(info))])
  set(this.gui.getHandle('sliderTrack'), 'Enable', 'On');
  
  sliderStep = [1 10] / (numel(info) - 1);
  set(this.gui.getHandle('sliderTrack'), 'Min', 1, 'Max', numel(info), 'SliderStep', sliderStep, 'Value', 1)
  set(this.gui.getHandle('buttonTrack'), 'Enable', 'on')
    
  
  %% set status to ready
  this.gui.changeStatus('statusMain', 'Ready...');


end

