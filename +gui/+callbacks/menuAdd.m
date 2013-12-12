function menuAdd(this)

  this.changeStatus('statusMain', 'Busy...');
  
  % read in filename
  [file, path] = uigetfile('*.tif*','Waehlen Sie eine Bilddatei',...
                                'MultiSelect','on');
    
  % check if uigetfile was cancled
  if path ~= 0 
    
    data.locate = [];
    data.track = [];
    data.calibrate = [];
    data.config = [];
    
    hTree = this.getHandle('treeMain');
    hTree.add(path, false, [], 'Workspace');
    hTree.add(file, true, data, path);
    
  end

  this.changeStatus('statusMain', 'Ready...');

end

