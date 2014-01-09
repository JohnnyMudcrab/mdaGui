function calibrate(this)
  
  hTree = this.gui.getHandle('treeMain');
  [selectedNodes, ~, ~] = hTree.getSelectedNodes();
  
  if isempty(selectedNodes)
    msgbox('No Image selected','No Image','Help')
    return
  end
  
  this.locate(false);
  this.track(false);

  voxelSize   = [this.gui.getText('textLocatePixelSizeX', 'numeric'),...
               this.gui.getText('textLocatePixelSizeY', 'numeric'),1];     
  fwhm_offset = this.gui.getText('textCalibrateFWHM', 'numeric');
  lens_shift  = this.gui.getText('textCalibrateLens', 'numeric');
  
  individualCalibration = this.gui.getValue('checkboxCalibrate');
  
  FWHMmin = this.gui.getText('textLocateFWHMmin', 'numeric');
  FWHMmax = this.gui.getText('textLocateFWHMmax', 'numeric');
  
  fwhmConfX = [0    0   0   0    0   -fwhm_offset  lens_shift];
  fwhmConfY = [0    0   0   0    0       0              0];

  
  % loop over nodes
  n = size(selectedNodes,1);

  for i = 1:1:n
    
    if(this.stop)
      this.stop = false;
      return;
    end
    
    this.gui.changeStatus('statusMain', ['Busy...Calibrating Image ' num2str(i) ' of ' num2str(n)]);

    data = selectedNodes{i}.handle.UserData;

    if ~(individualCalibration)

      % combine tracks
      temp  = data.track.trajectories(:,1);
      combined_p = {};
      
      for j = 1:size(temp,1)
        combined_p = [combined_p; temp{j,1}(:,1:2); {[] []}]; %#ok<AGROW>
      end

      clear temp

      if isempty(combined_p)
        fwhmConfX = [];
        fwhmConfY = [];
        success = -1;
        control = {};
      else
        [~, fwhmConfX, fwhmConfY, success, control] = ...
          bin.calibrateInvoker(combined_p, FWHMmin, FWHMmax, fwhmConfX, fwhmConfY, lens_shift, voxelSize, 1);

      end
      
      
    end

    m = size(data.track.trajectories,1);
    p = cell(n,1);

    for j = 1:1:m
      
      this.gui.changeStatus('statusMain', ['Busy...Calibrating Image ' num2str(i) ' of ' num2str(n) ...
                            ', Track ' num2str(j) ' of ' num2str(m)]);

      if(this.stop)
        this.stop = false;
        return;
      end

      if(~individualCalibration && success >= 0)
        [p_current, fwhmConfX, fwhmConfY, ~, control] = ...
          bin.calibrateInvoker(data.track.trajectories{j}, FWHMmin, FWHMmax, fwhmConfX, fwhmConfY, lens_shift, voxelSize, 0);
      else
        if (~individualCalibration && success == -1)
          p_current = data.track.trajectories{j};
        end
      end

      if(individualCalibration)
        [p_current, fwhmConfX, fwhmConfY, ~, control] = ...
          bin.calibrateInvoker(data.track.trajectories{j}, FWHMmin, FWHMmax, fwhmConfX, fwhmConfY, lens_shift, voxelSize, 1);
      end

      p{j,1} = p_current;
      p{j,2} = fwhmConfX;
      p{j,3} = fwhmConfY;
      p{j,4} = data.track.trajectories{j,2};
      p{j,5} = data.track.trajectories{j,3};
      
      if j == 1
        p{j,6} = control;
      else
        p{j,6} = {};
      end

    end

    data.calibrate = p;
    selectedNodes{i}.handle.UserData = data;

  end

end

