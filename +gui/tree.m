classdef tree < handle
  
  properties(Access = private)
    
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
    
    % getter
    function selectedNodes = getSelectedNodes(this)
      selectedNodes = this.handle.getSelectedNodes;
    end
    
    % utilities
    function disableMultipleSelection(this)
      
      this.handle.setMultipleSelectionEnabled(false);
      
    end
    
    function enableMultipleSelection(this)
      
      this.handle.setMultipleSelectionEnabled(true);
      
    end
    
    function saveExpansionState(this)
      
    end
    
    function restoreExpansionState(this)
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

