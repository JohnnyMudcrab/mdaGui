classdef mda < handle

  properties (Access = public)
    
    gui;
    stop;
    
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
      
    end
    
    function locate(this)
      
      hTree = this.gui.getHandle('treeMain');
      
      
      [selectedNodes, filenames, ~] = hTree.getSelectedNodes();
      
      if isempty(selectedNodes)
        msgbox('No Image selected','No Image','Help')
        return
      end
      
      n = numel(selectedNodes);
      
      for i = 1:n
        
        filename = filenames{i}; 
        
        info = imfinfo(filename); 
        nFrames = numel(info); 

        currentNode = selectedNodes{i};

        data = currentNode.handle.UserData;

        data.locate = cell(nFrames,1); 


        % Locate Particles in every Frame

        for j = 1:1:nFrames 
          
          if(this.stop)
            this.stop = false;
            return;
          end
          
          this.gui.changeStatus('statusMain', ['Busy...Processing Image ' num2str(i) ' of ' num2str(n) ...
                                ', Frame ' num2str(j) ' of ' num2str(nFrames)]);
                              
          drawnow

          img = imread(filename,'Index',j,'Info',info); 

          gConfs = locateParticles2D(img, []); 

          data.locate{j,1} = gConfs;

        end

        currentNode.handle.UserData = data;
        
      end
      
    end
    
  end
  
  methods(Static)

    function closingFcn(~,~)
      
      hMda = getappdata(0, 'hMda');
      
      if ~isempty(hMda)
        delete(hMda.gui)
        delete(hMda)
        rmappdata(0, 'hMda');
      end

      delete(gcf);
      
    end
  
  end
  
end

