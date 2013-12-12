function buttonInfo(this)

  %set status to busy
  this.changeStatus('statusMain', 'Busy...');
  
  %get tree
  hTree = this.getHandle('treeMain');
  
  selectedNodes = hTree.getSelectedNodes();
  
  % non selected
  if isempty(selectedNodes)
    msgbox('No Image selected','No Image','Help')
    return
  end
  
  % more then one selected
  if numel(selectedNodes) > 1 || ~selectedNodes(1).isLeaf
    msgbox('Too many Images selected','Too many Images','Help')
    return
  end
  
  % get file location
  file = selectedNodes(1).handle.UserData.string;
  path = selectedNodes(1).getParent.handle.UserData.string;
  
  % show picture in axes
  axes(this.getHandle('axesInfo'))
  info = imfinfo([path file]); 
  img = imread([path file],'Index',1,'Info',info); 
  imshow(double(img) / 255) 
    
  % update info panel
  set(this.getHandle('textInfoName2'), 'String', [path file])
  set(this.getHandle('textInfoResolution2'), 'String', [num2str(info(1, 1).Width) ' x ' num2str(info(1, 1).Height)])
  set(this.getHandle('textInfoNumber2'), 'String', num2str(numel(info)))
  set(this.getHandle('textInfoFrame'), 'String', ['1 / ' num2str(numel(info))])
  
  % enable slider
  set(this.getHandle('sliderInfo'), 'Enable', 'On');
  sliderStep = [1 10] / (numel(info) - 1);
  set(this.getHandle('sliderInfo'), 'Min', 1, 'Max', numel(info), 'SliderStep', sliderStep, 'Value', 1)
    
  %set status to ready
  this.changeStatus('statusMain', 'Ready...');

end

