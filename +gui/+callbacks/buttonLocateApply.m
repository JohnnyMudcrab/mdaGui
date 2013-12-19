function buttonLocateApply(this)

  this.changeStatus('statusMain', 'Busy...');
  
  path = this.getText('textInfoName2');
  index = int32(get(this.getHandle('sliderLocate'), 'Value'));
  
  img = imread(path,'Index',index);
  
  gConfs = locateParticles2D(img, []); 
  
  % gConfs(7) = x; gConfs(8) = y
  axes(this.getHandle('axesLocate'))
  
  imshow(double(img) / 255)
  hold on
  
  for i = 1:size(gConfs,1)
    rectangle('Position',[gConfs(i,7) - 10,gConfs(i,8) - 10,20,20],...
              'Curvature',[1,1],'EdgeColor','r')
  end

  hold off
  
  this.changeStatus('statusMain', 'Ready...');
  
end
