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
      
      % disable button
      set(this.gui.getHandle('buttonPreview'), 'Enable', 'Off')
      set(this.gui.getHandle('buttonLocate'), 'Enable', 'Off')
      set(this.gui.getHandle('buttonCreateMask'), 'Enable', 'Off')
      set(this.gui.getHandle('buttonDeleteMask'), 'Enable', 'Off')
      set(this.gui.getHandle('buttonTrack'), 'Enable', 'Off')
      set(this.gui.getHandle('buttonCalibrate'), 'Enable', 'Off')
      set(this.gui.getHandle('menuStop'), 'Enable', 'Off')
      
      % text format
      set(this.gui.getHandle('textLocateState'), 'FontWeight', 'bold', 'ForegroundColor', [1 0 0])
      set(this.gui.getHandle('textTrackState'), 'FontWeight', 'bold', 'ForegroundColor', [1 0 0])
      set(this.gui.getHandle('textCalibrateState'), 'FontWeight', 'bold', 'ForegroundColor', [1 0 0])
      
      % set mouseclick action for hTree
      hTree = this.gui.getHandle('treeMain');
      jTree = handle(hTree.handle.getTree,'CallbackProperties');
      set(jTree, 'MousePressedCallback', @this.mousePressedCallback);    
      
      % listbox
      hListbox = this.gui.getHandle('listCalibrate');
      set(hListbox, 'Min', 1, 'Max', 5);
      set(hListbox, 'Callback', @(src, event)gui.callbacks.listboxCalibrate(this));
      
      
    end
    
    function update(this)
      
      data = this.currentNode.handle.UserData;
      
      % update state
      if ~isempty(data.locate)
        set(this.gui.getHandle('textInfoStateLocate'), 'String', 'Located')
        set(this.gui.getHandle('textLocateState'), 'String', 'Located', 'ForegroundColor', [0 1 0])
      else
        set(this.gui.getHandle('textInfoStateLocate'), 'String', 'Not Located Yet')
        set(this.gui.getHandle('textLocateState'), 'String', 'Not Located', 'ForegroundColor', [1 0 0])
      end

      if ~isempty(data.track)
        set(this.gui.getHandle('textInfoStateTrack'), 'String', 'Tracked')
        set(this.gui.getHandle('textLocateTrack'), 'String', 'Tracked', 'ForegroundColor', [0 1 0])
      else
        set(this.gui.getHandle('textInfoStateTrack'), 'String', 'Not Tracked Yet')
        set(this.gui.getHandle('textTrackState'), 'String', 'Not Tracked', 'ForegroundColor', [1 0 0])
      end
      
      if ~isempty(data.calibrate)
        set(this.gui.getHandle('textInfoStateCalibrate'), 'String', 'Calibrated')
        set(this.gui.getHandle('textCalibrateState'), 'String', 'Calibrated', 'ForegroundColor', [0 1 0])
      else
        set(this.gui.getHandle('textInfoStateCalibrate'), 'String', 'Not Calibrated Yet')
        set(this.gui.getHandle('textCalibrateState'), 'String', 'Not Calibrated', 'ForegroundColor', [1 0 0])
      end
      
      % update calibrate list
      hListbox = this.gui.getHandle('listCalibrate');
      data = sortrows(data.calibrate,-5);
      
      n = size(data,1);
      
      list = cell(n,1);
      
      for i = 1:n
        
        list{i} = ['Track ID ' num2str(data{i,4}) ' (' num2str(data{i,5}) ')'];
        
      end
      
      set(hListbox, 'String', list);
    
    end

    
    apply(this)
    
    calibrate(this)
    
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

