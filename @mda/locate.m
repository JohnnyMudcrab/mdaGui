function locate(this, overwrite)
      
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
    
    if isempty(data.locate) || overwrite

      data.locate = cell(nFrames,1); 


      % Locate Particles in every Frame

      for j = 1:1:nFrames 

        if(this.stop)
          this.stop = false;
          return;
        end

        this.gui.changeStatus('statusMain', ['Busy...Locating Image ' num2str(i) ' of ' num2str(n) ...
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
