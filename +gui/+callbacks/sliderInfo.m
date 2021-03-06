function sliderInfo(this)

  % get image location
  path = get(this.getHandle('textInfoName2'), 'String');
  
  % get slider value
  index = int32(get(this.getHandle('sliderInfo'),'Value'));
  set(this.getHandle('sliderInfo'),'Value',index);
  set(this.getHandle('textInfoFrame'), 'String', ...
    [num2str(index) ' / ' get(this.getHandle('textInfoNumber2'), 'String')])
  
  % show picture in axes
  axes(this.getHandle('axesInfo'))
  set(this.getHandle('axesInfo'), 'Nextplot','replacechildren')
  img = imread(path,'Index',index); 
  imshow(img, [this.getText('textInfoHistMin', 'numeric') this.getText('textInfoHistMax', 'numeric')]) 
  
  hMda = getappdata(0,'hMda');
  data = hMda.currentNode.handle.UserData;
  
  % draw ROI
  if(isfield(data,'mask'))
    hMda.drawROI(data.mask.position)
  end
  
  
end

