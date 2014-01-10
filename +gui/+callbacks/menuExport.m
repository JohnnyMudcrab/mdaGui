function menuExport(this)

  this.changeStatus('statusMain', 'Busy...');
  
  if isunix
    
    hTree = this.getHandle('treeMain');


    [selectedNodes, ~, ~] = hTree.getSelectedNodes();

    if isempty(selectedNodes)
      msgbox('No Image selected','No Image','Help')
      return
    end

    n = numel(selectedNodes);
    
    data = cell(n,1);
    
    for i = 1:n
      
      data{i} = selectedNodes{i}.handle.UserData;
      
    end
    
    assignin('base','data', data)
    
  else
    
    hTree = this.getHandle('treeMain');


    [selectedNodes, ~, ~] = hTree.getSelectedNodes();

    if isempty(selectedNodes)
      msgbox('No Image selected','No Image','Help')
      return
    end
    
    if numel(selectedNodes) > 1
      msgbox('Please only select one movie. The export of more then one movie is currently not supported!', ...
             'Sorry', 'Help')
      return
    end
    
    data = selectedNodes{1}.handle.UserData.calibrate;
    
    if isempty(data)
      msgbox('Please calibrate befor export', '', 'Help')
      return
    else
    
    n = size(data,1);
    m = 0;
    
    for i = 1:n
      
      this.changeStatus('statusMain', ['Busy ... Calculation Size of Excel Sheet for Track ' num2str(i) ' of ' num2str(n)]);
      
      if size(data{i,1}{1,1},2) > 9
        m = m + size(data{i,1},1);
      end
      
    end
    
    sheet = cell(m, 5);
    
    iterator = 1;
     
    for i = 1:n
      
      m = size(data{i,1},1);
      
      this.changeStatus('statusMain', ['Busy ... Exporting Track ' num2str(i) ' of ' num2str(n)]);
      
      for j = 1:m
        
        sheet{iterator, 1} = data{i, 4};              % track
        sheet{iterator, 2} = data{i, 1}{j,2};         % frame
        sheet{iterator, 3} = data{i, 1}{j,1}(1,7);    % x
        sheet{iterator, 4} = data{i, 1}{j,1}(1,8);    % y
        sheet{iterator, 5} = data{i, 1}{j,1}(1,11);   % z
        
        iterator = iterator + 1;
        
      end
      
    end
    
    [file, path] = uiputfile('*.xls');
    
      
    xlswrite([path file], sheet);
      
    end
        
  end
  
  this.changeStatus('statusMain', 'Ready...');

end

