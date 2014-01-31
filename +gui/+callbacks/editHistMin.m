function editHistMin(this)

  path = get(this.gui.getHandle('textInfoName2'), 'String');

  % show picture in info
  axes(this.gui.getHandle('axesInfo'))
  index = int32(get(this.gui.getHandle('sliderInfo'),'Value'));
  img = imread(path,'Index',index); 
  set(this.gui.getHandle('axesInfo'), 'Nextplot','replacechildren')
  imshow(img, [this.gui.getText('textInfoHistMin', 'numeric') this.gui.getText('textInfoHistMax', 'numeric')])
  
  % show picture in locate
  axes(this.gui.getHandle('axesLocate'))
  index = int32(get(this.gui.getHandle('sliderLocate'),'Value'));
  img = imread(path,'Index',index); 
  set(this.gui.getHandle('axesLocate'), 'Nextplot','replacechildren')
  imshow(img, [this.gui.getText('textInfoHistMin', 'numeric') this.gui.getText('textInfoHistMax', 'numeric')]) 
  
  % show picture in locate
  axes(this.gui.getHandle('axesTrack'))
  index = int32(get(this.gui.getHandle('sliderTrack'),'Value'));
  img = imread(path,'Index',index); 
  set(this.gui.getHandle('axesTrack'), 'Nextplot','replacechildren')
  imshow(img, [this.gui.getText('textInfoHistMin', 'numeric') this.gui.getText('textInfoHistMax', 'numeric')]) 


end