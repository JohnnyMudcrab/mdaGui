function sliderLocate(this)

  hMda = getappdata(0,'hMda');
  
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
      rectangle('Position',[data.locate{index}(i,7) - 10,data.locate{index}(i,8) - 10,20,20],...
                'Curvature',[1,1],'EdgeColor','r')
    end
    hold off
    
  end
  
end

