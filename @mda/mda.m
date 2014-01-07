classdef mda < handle

  properties (Access = public)
    
    gui;
    stop;
    currentNode;
    
  end
  
  methods (Access = public)
    
    function this = mda()
      
      this.stop = false;
       
    end
    
    function init(this)
      
      this.gui = gui.framework('+gui/gui.cfg','gui'); %#ok<*PROP>
      
      % disable slider
      set(this.gui.getHandle('sliderInfo'), 'Enable', 'Off')
      set(this.gui.getHandle('sliderLocate'), 'Enable', 'Off')
      set(this.gui.getHandle('sliderTrack'), 'Enable', 'Off')
      set(this.gui.getHandle('sliderCalibrate'), 'Enable', 'Off')
      
      % disable button
      set(this.gui.getHandle('buttonLocateApply'), 'Enable', 'Off')
      set(this.gui.getHandle('buttonLocateApplyAll'), 'Enable', 'Off')
      set(this.gui.getHandle('buttonTrack'), 'Enable', 'Off')
      
      % set mouseclick action for hTree
      hTree = this.gui.getHandle('treeMain');
      jTree = handle(hTree.handle.getTree,'CallbackProperties');
      set(jTree, 'MousePressedCallback', @this.mousePressedCallback);    
      
    end

    
    apply(this)
    
    locate(this, overwrite)
    
    track(this, overwrite)
    
  end
  
  methods(Static)
    
    function mousePressedCallback(~, eventData)
      
      this = getappdata(0, 'hMda');

      if eventData.getClickCount==2 
        this.apply()
      end
      
    end     


    function closingFcn(~,~)
      
      hMda = getappdata(0, 'hMda');
      
      if ~isempty(hMda)
        delete(hMda.gui)
        delete(hMda)
        rmappdata(0, 'hMda');
      end

      delete(gcf);
      
    end
    
    function[] = drawROI(position)

      hold on
      n = size(position,1);
      for i = 1:n-1
        line([position(i,1), position(i+1,1)],[position(i,2), position(i+1,2)],'Color',[1 0 0])
      end
      line([position(1,1), position(n,1)],[position(1,2), position(n,2)],'Color',[1 0 0]) 
      hold off

    end
    

  
  end
  
end

