function menuConvert(this)
  [file,path,~] = uigetfile('*.mat','Select the MATLAB code file');
  
  if path ~=0
    
    filename = textscan(file, '%s %s', 'delimiter','.'); 

    % check if file is a mat-file
    if strcmp(filename{2}{1},'mat') == 0
       waitfor(msgbox('Kein mat ausgewaehlt!',...
                   'ERROR','error')); 
       return; 
    end 

    % load file
    load([path file]);
  
    n = size(saveList,1); 
    
    for i = 1:1:n
      
      result = strfind(saveList{i,1},'/');
      file = saveList{i,1}(result(end) + 1: size(saveList{i,1},2));
      path = saveList{i,1}(1:result(end));
      
      data.locate = saveList{i,2}.locate;
      data.track = saveList{i,2}.track;
      data.calibrate = saveList{i,2}.calibrate;
      data.mask = saveList{i,2}.mask;
      data.string = [path file];
      
      
      saveList{i,2} = data;
      
    end
    
    [file, path] = uiputfile;
    save([path file], 'saveList','treeKey')
  
    
  end
  
  
  
end

