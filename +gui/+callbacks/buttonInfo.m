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
  
  % update info panel
  set(this.getHandle('textInfoName2'), 'String', [path file])
  
  % show picture in axes
  axes(this.getHandle('axesInfo'))
  info = imfinfo([path file]); 
  img = imread([path file],'Index',1,'Info',info); 
  imshow(double(img) / 255) 
  
  % enable slider
  set(this.getHandle('sliderInfo'), 'Enable', 'On');
    
  %set status to ready
  this.changeStatus('statusMain', 'Ready...');

end

