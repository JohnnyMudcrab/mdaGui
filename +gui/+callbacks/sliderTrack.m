function sliderTrack(this)

  hMda = getappdata(0,'hMda');
  
  % get image location
  path = get(this.getHandle('textInfoName2'), 'String');
  
  % get slider value
  index = int32(get(this.getHandle('sliderTrack'),'Value'));
  set(this.getHandle('sliderTrack'),'Value',index);
  set(this.getHandle('textTrackFrame'), 'String', ...
    [num2str(index) ' / ' get(this.getHandle('textInfoNumber2'), 'String')])
  
  % show picture in axes
  axes(this.getHandle('axesTrack'))
  img = imread(path,'Index',index); 
  imshow(double(img) / 255) 
  
  data = hMda.currentNode.handle.UserData;
  
  % draw ROI
  if(isfield(data,'mask'))
    hMda.drawROI(data.mask.position)
  end
  
  % show results
  if ~isempty(data.locate)
    
    hold on
    for i = 1:size(data.locate{index},1)
      rectangle('Position',[data.locate{index}(i,7) - 10,data.locate{index}(i,8) - 10,20,20],...
                'Curvature',[1,1],'EdgeColor','r')
              
      if ~isempty(data.track)
        text(data.track.frames{index}(i,8),data.track.frames{index}(i,9) - 5,num2str(data.track.frames{index}(i,1)))
      end
      
    end
    hold off
    
  end
  
end

