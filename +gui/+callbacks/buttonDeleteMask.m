function buttonDeleteMask(this)

  hMda = getappdata(0,'hMda');
  
  data = hMda.currentNode.handle.UserData;
  
  if(isfield(data,'mask'))
    data = rmfield(data,'mask');
  end
  
  hMda.currentNode.handle.UserData = data;
  
end

