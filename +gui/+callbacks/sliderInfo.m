function sliderInfo(this)

  % get image location
  path = get(this.getHandle('textInfoName2'), 'String');
  
  % get slider value
  index = int32(get(this.getHandle('sliderInfo'),'Value'));
  set(this.getHandle('sliderInfo'),'Value',index);
  set(this.getHandle('textInfoFrame'), 'String',[num2str(index) ' / '])
  
  % show picture in axes
  axes(this.getHandle('axesInfo'))
  img = imread(path,'Index',index); 
  imshow(double(img) / 255) 
  
end
