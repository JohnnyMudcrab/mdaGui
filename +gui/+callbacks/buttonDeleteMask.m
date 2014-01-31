function buttonDeleteMask(this)

  hMda = getappdata(0,'hMda');
  
  data = hMda.currentNode.handle.UserData;
  
  if(isfield(data,'mask'))
    data = rmfield(data,'mask');
  end
  
  hMda.currentNode.handle.UserData = data;
  
  % get path of image
  path = get(this.getHandle('textInfoName2'), 'String');
  
  % show picture in locate
  axes(this.getHandle('axesLocate'))
  index = int32(get(this.getHandle('sliderLocate'),'Value'));
  img = imread(path,'Index',index); 
  set(this.getHandle('axesLocate'), 'Nextplot','replacechildren')
  imshow(img, [this.getText('textInfoHistMin', 'numeric') this.getText('textInfoHistMax', 'numeric')]) 
  
  % show picture in locate
  axes(this.getHandle('axesTrack'))
  index = int32(get(this.getHandle('sliderTrack'),'Value'));
  img = imread(path,'Index',index); 
  set(this.getHandle('axesTrack'), 'Nextplot','replacechildren')
  imshow(img, [this.getText('textInfoHistMin', 'numeric') this.getText('textInfoHistMax', 'numeric')]) 
  
end

