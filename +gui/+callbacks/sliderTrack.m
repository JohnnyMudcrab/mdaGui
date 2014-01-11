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
  set(this.getHandle('axesTrack'), 'Nextplot','replacechildren')
  imshow(img, [this.getText('textInfoHistMin', 'numeric') this.getText('textInfoHistMax', 'numeric')]) 
  
  data = hMda.currentNode.handle.UserData;
  
  % draw ROI
  if(isfield(data,'mask'))
    hMda.drawROI(data.mask.position)
  end
  
  % show results
  if ~isempty(data.locate)
    
    hold on
    for i = 1:size(data.locate{index},1)
%       rectangle('Position',[data.locate{index}(i,7) - 10,data.locate{index}(i,8) - 10,20,20],...
%                 'Curvature',[1,1],'EdgeColor','r')
              
      drawGaussian([data.locate{index}(i,7), data.locate{index}(i,8), ...
              2 * data.locate{index}(i,3) / 2.3548200450309, ...
              2 * data.locate{index}(i,4) / 2.3548200450309, ...
              0], 'g');
              
      if ~isempty(data.track)
        text(data.track.frames{index}(i,8),data.track.frames{index}(i,9) - 5, ...
             num2str(data.track.frames{index}(i,1)), 'Color', [1 0 0])
      end
      
    end
    hold off
    
  end
  
end

