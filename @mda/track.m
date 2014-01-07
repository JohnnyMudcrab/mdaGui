function track(this, overwrite)

  hTree = this.gui.getHandle('treeMain');


  [selectedNodes, ~, ~] = hTree.getSelectedNodes();

  if isempty(selectedNodes)
    msgbox('No Image selected','No Image','Help')
    return
  end
  
  this.locate(false);


  %% Loop over unique Nodes
  
  n = numel(selectedNodes);

  for i = 1:1:n
    
    this.gui.changeStatus('statusMain', ['Busy...Tracking Image ' num2str(i) ' of ' num2str(n)]);
    drawnow

    data = selectedNodes{i}.handle.UserData;


    % Connect Particles
    gConfs = data.locate;
    m = size(gConfs,1);

    for j = 1:1:m
      if ~isempty(gConfs{j})
        gConfs{j}(:,1:6) = [];
      end
    end
    
    maxDist = this.gui.getText('textTrackMaxDisp2', 'numeric');
    maxLost = this.gui.getText('textTrackMaxLost2', 'numeric');

    [tracks, adjacency_tracks] = simpletracker(gConfs,...
      'MaxLinkingDistance', maxDist, ...
      'MaxGapClosing', maxLost);


    n_tracks = numel(tracks);

    frames = cell(m,1);
    trajectories = cell(n_tracks,1);
    iterator = 1;

    for j = 1:1:n_tracks

        o = size(tracks{j},1);

        trajectories{j,2} = j;
        trajectories{j,3} =  numel(tracks{j}) - sum(isnan(tracks{j}));

        for k = 1:1:o
          
          if ~isnan(tracks{j}(k))
            frames{k,1}(size(frames{k,1},1)+1,1) = j;
            frames{k,1}(size(frames{k,1},1),2:9) = data.locate{k,1}(tracks{j}(k),1:8);

            trajectories{j,1}{iterator,1} = data.locate{k,1}(tracks{j}(k),1:8);
            trajectories{j,1}{iterator,2} = k;
            iterator = iterator + 1;
          end
          
        end

      iterator = 1;

    end

    data.track = [];

    data.track.frames = frames;
    data.track.trajectories = trajectories;
    
    selectedNodes{i}.handle.UserData = data;    

  end 




end

