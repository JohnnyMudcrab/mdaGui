function menuAdd(this)

  this.changeStatus('statusMain', 'Busy...');
  
  % read in filename
  [file, path] = uigetfile('*.tif*','Waehlen Sie eine Bilddatei',...
                                'MultiSelect','on');
    
  % check if uigetfile was cancled
  if path ~= 0 
    
    hTree = this.getHandle('treeMain');
    hTree.add(path, false, [], 'Workspace');
    hTree.add(file, true, [], path);
    
  end

  this.changeStatus('statusMain', 'Ready...');

end

