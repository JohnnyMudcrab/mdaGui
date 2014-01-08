function buttonPreview(this)

  this.changeStatus('statusMain', 'Busy...');
  
  % disable gui while busy
  this.disableGui();
  set(this.getHandle('menuStop'), 'enable', 'on')
  
  % locate particles in current picture
  path = this.getText('textInfoName2');
  index = int32(get(this.getHandle('sliderLocate'), 'Value'));
  
  img = imread(path,'Index',index);
  
  hMda = getappdata(0, 'hMda');
  data = hMda.currentNode.handle.UserData;
  
  if(isfield(data,'mask'))
    gConfs = locateParticles2D(img, data.mask.bw); 
  else
    gConfs = locateParticles2D(img, []); 
  end
  
  % gConfs(7) = x; gConfs(8) = y
  axes(this.getHandle('axesLocate'))
  
  imshow(double(img) / 255)
  hold on
  
  for i = 1:size(gConfs,1)
    rectangle('Position',[gConfs(i,7) - 10,gConfs(i,8) - 10,20,20],...
              'Curvature',[1,1],'EdgeColor','r')
  end

  hold off

  if(isfield(data,'mask'))
    hMda.drawROI(data.mask.position)
  end
  
  
  % enable gui again
  this.enableGui();
  
  this.changeStatus('statusMain', 'Ready...');
  
end

