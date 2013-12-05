classdef mda < handle

  properties (Access = public)
    
    gui;
    
  end
  
  methods (Access = public)
    
    function this = mda()
       
    end
    
    function init(this)
      
      this.gui = gui.framework('+gui/gui.cfg','gui'); %#ok<*PROP>
      
      set(this.gui.getHandle('sliderInfo'), 'Enable', 'Off')
      
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

