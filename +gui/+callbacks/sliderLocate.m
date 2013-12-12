function sliderLocate(this)
  
% get image location
  path = get(this.getHandle('textInfoName2'), 'String');
  
  % get slider value
  index = int32(get(this.getHandle('sliderLocate'),'Value'));
  set(this.getHandle('sliderLocate'),'Value',index);
  set(this.getHandle('textLocateFrame'), 'String', ...
    [num2str(index) ' / ' get(this.getHandle('textInfoNumber2'), 'String')])
  
  % show picture in axes
  axes(this.getHandle('axesLocate'))
  img = imread(path,'Index',index); 
  imshow(double(img) / 255) 
  
end

