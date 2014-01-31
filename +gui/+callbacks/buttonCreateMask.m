function buttonCreateMask(this)

  hMda = getappdata(0,'hMda');

  h = impoly(this.getHandle('axesLocate'));
  position = wait(h);
  bw = createMask(h);
  
  data = hMda.currentNode.handle.UserData;
  data.mask.bw = bw;
  data.mask.position = position;
  hMda.currentNode.handle.UserData = data;
  
  delete(h);
  

  % show picture in locate
  axes(this.getHandle('axesLocate'))
  hMda.drawROI(position);
  
  % show picture in locate
  axes(this.getHandle('axesTrack'))
  hMda.drawROI(position);
  
  
end

