classdef tree < handle
  
  properties(Access = public)
    
    handle
    root
    container
    
  end
  
  methods(Access = public)
    
    function this = tree(string, parent)
      
      this.root = uitreenode('v0',string,string,[], false); 
      [this.handle, this.container] = uitree('v0', 'Root',this.root);
      
      set(this.container, 'Parent', parent);
      
    end
    
    
    % add and remove nodes
    function add(this, nodes, isLeaf, data, parent)
      
      if ~iscell(nodes)
        nodes = cellstr(nodes);
      end
      
      parentNode = this.findParent(parent);
      
      for i = 1:numel(nodes)
        
        if ~this.isDuplicate(nodes{i}, parentNode)
          data.string = nodes{i};
          if numel(data.string) > 20
            node = uitreenode('v0',handle(this.handle), ['...' data.string(end-20:end)], [], isLeaf); %#ok<*PROP>
          else
            node = uitreenode('v0',handle(this.handle), data.string, [], isLeaf); %#ok<*PROP>
          end
          
          node.handle.UserData = data;
          parentNode.add(node);
        end
        
      end
      
      this.handle.reloadNode(this.root);
       
    end
    
    function removeSelectedNodes(this)
      
      selectedNodes = this.getSelectedNodes(); 

      if ~isempty(selectedNodes)

        for i = 1:1:size(selectedNodes,1)

          level = selectedNodes{1}.getLevel; 

          if level 

            parent = selectedNodes{1}.getParent;
            selectedNodes{1}.removeFromParent; 
            if ~parent.getChildCount && ~isempty(parent.getParent)
              this.root.remove(parent);
            end      

          else 
            this.root.removeAllChildren; 
          end  
          
        end
      end

      this.reloadTree(); 
      
    end
    
    % getter
    function [selectedNodes, filenames, error] = getSelectedNodes(this, option)
      
      if nargin < 2
        option = 'multiple';
      end
      
      error = [];
      
      selectedNodes = this.handle.getSelectedNodes;
      
      
      % check if empty
      if isempty(selectedNodes)
        error = 1;
        return
      end
      
      
      % check if only one leaf is selected
      if strcmp(option, 'single')
        
        if numel(selectedNodes) > 1 || ~selectedNodes(1).isLeaf
          selectedNodes = [];
          error = 2;
        end
        
        return
        
      end

      
      % make it unique    
      n = size(selectedNodes,1);
      
      filenames  = {};
      nodes   = {};
      count   = 0;

      for i = 1:1:n
        level = selectedNodes(1).getLevel;
        if level == 0
          tempParent = selectedNodes(1).getFirstChild;
          while(~isempty(tempParent))
            tempChild = tempParent.getFirstChild;
            while(~isempty(tempChild))
              count = count + 1;
              nodes{count,1} = tempChild;
              filenames{count,1} = [char(tempParent.handle.UserData.string) char(tempChild.handle.UserData.string)];
              tempChild = tempParent.getChildAfter(tempChild);
            end
            tempParent = selectedNodes(1).getChildAfter(tempParent);
          end
        elseif level == 1
          tempChild = selectedNodes(1).getFirstChild;
          while(~isempty(tempChild))
            tempParent = tempChild.getParent;
            count = count + 1;
            nodes{count,1} = tempChild;
            filenames{count,1} = [char(tempParent.handle.UserData.string) char(tempChild.handle.UserData.string)];
            tempChild = selectedNodes(1).getChildAfter(tempChild);
          end
        else
           tempChild = selectedNodes(1);
           tempParent = tempChild.getParent;
           count = count + 1;
           nodes{count,1} = tempChild;
           filenames{count,1} = [char(tempParent.handle.UserData.string) char(tempChild.handle.UserData.string)];
        end
      end

     [filenames, nodeIndex] = unique(filenames,'first');

     selectedNodes = nodes(nodeIndex);

 
    end
    
    % utilities
    function disableMultipleSelection(this)
      
      this.handle.setMultipleSelectionEnabled(false);
      
    end
    
    function enableMultipleSelection(this)
      
      this.handle.setMultipleSelectionEnabled(true);
      
    end
    
    function save(this)
      
      [file, path] = uiputfile;
      
      if path ~= 0 
        
        count = this.root.LeafCount;
        
        if count
          
          saveList = cell(count, 2);
          listCount = 0;
          tempParent = this.root.getFirstChild;
          
          while(~isempty(tempParent))
            tempChild = tempParent.getFirstChild;
            
            while(~isempty(tempChild))
              listCount = listCount + 1;
              saveList{listCount,1} = [char(tempParent.handle.UserData.string) char(tempChild.handle.UserData.string)];
              saveList{listCount,2} = tempChild.handle.UserData; 
              tempChild = tempParent.getChildAfter(tempChild);
            end
            
            tempParent = this.root.getChildAfter(tempParent);
          end
          
        else
          saveList = {}; %#ok<NASGU>
        end

        treeKey = 1; %#ok<NASGU>
        save([path file], 'saveList','treeKey')
        
      end
      
    end
    
    function load(this)
      
        [file,path,~] = uigetfile('*.mat','Select the MATLAB code file');

        if path ~= 0 

          filename = textscan(file, '%s %s', 'delimiter','.'); 
                                                             
          % check if file is a mat-file
          if strcmp(filename{2}{1},'mat') == 0
             waitfor(msgbox('Kein mat ausgewaehlt!',...
                         'ERROR','error')); 
             return; 
          end

          % load file
          load([path file]);
        end
        
        % check if it is a tree file
        if ismember('treeKey',who) 

          this.root.removeAllChildren;

          n = size(saveList,1); %#ok<USENS>

          % loop over tree
          for i = 1:1:n
            result = strfind(saveList{i,1},'/');
            file = saveList{i,1}(result(end) + 1: size(saveList{i,1},2));
            path = saveList{i,1}(1:result(end));
            userData = saveList{i,2};

            this.add(path, false, [], 'Workspace');
            this.add(file, true, userData, path);
            
          end

          this.reloadTree();
        else
           warndlg('Not a tree file');
        end

      
    end
    
    function saveExpansionState(this)
      
    end
    
    function restoreExpansionState(this)
    end
    
    function reloadTree(this)
      this.handle.reloadNode(this.root); 
    end
    
  end
  
  
  methods(Access = private)
    
    function parentNode = findParent(this, parent)
      
      parentNode = [];
      
      if strcmp(this.root.getName, parent)
        parentNode = this.root;
      else
        
        if this.root.getChildCount
          tempNode = this.root.getFirstChild;
        else
          tempNode = [];
        end
        
        while(~isempty(tempNode))
          if strcmp(tempNode.handle.UserData.string,parent)
            parentNode = tempNode;
            tempNode = [];
          else
            tempNode = this.root.getChildAfter(tempNode);
          end
        end
        
      end
      
    end
     
  end
  
  
  methods(Static)

    function dublicate = isDuplicate(node, parent)
      
      dublicate = false;
      
      if parent.getChildCount
        tempNode = parent.getFirstChild;
      else
        tempNode = [];
      end

      while(~isempty(tempNode))
        if strcmp(tempNode.handle.UserData.string, node)
          dublicate = true;
          return
        else
          tempNode = parent.getChildAfter(tempNode);
        end
      end
      
    end
    
  end
  
end

