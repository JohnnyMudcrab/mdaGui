function sliderInfo(this)

  % get image location
  path = get(this.getHandle('textInfoName2'), 'String');
  
  % show picture in axes
  axes(this.getHandle('axesInfo'))
  img = imread(path,'Index',get(this.getHandle('sliderInfo'),'Value')); 
  imshow(double(img) / 255) 
  
end

