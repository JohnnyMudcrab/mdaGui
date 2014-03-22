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
      mask = data.mask.bw;
  else
      mask = [];
  end

  gConfs = bin.locateInvoker(img, ...
                         this.getText('textLocateTreshold2', 'numeric'), ...
                         this.getText('textLocatePixelSizeX', 'numeric'), ...
                         this.getText('textLocatePixelSizeY', 'numeric'), ...
                         this.getText('textLocateFWHMmin', 'numeric'), ...
                         this.getText('textLocateFWHMmax', 'numeric'), ...
                         mask);
  
  % gConfs(7) = x; gConfs(8) = y
  axes(this.getHandle('axesLocate'))
  
  imshow(img, [this.getText('textInfoHistMin', 'numeric') this.getText('textInfoHistMax', 'numeric')])
  hold on
  
  for i = 1:size(gConfs,1)            
    drawGaussian([gConfs(i,7), gConfs(i,8), ...
        2 * gConfs(i,3) / 2.3548200450309, ...
        2 * gConfs(i,4) / 2.3548200450309, ...
        0], 'g');
  end

  hold off

  if(isfield(data,'mask'))
    hMda.drawROI(data.mask.position)
  end
  
  
  % enable gui again
  this.enableGui();
  
  this.changeStatus('statusMain', 'Ready...');
  
end

