function listboxCalibrate(this)

      
      type = get(this.gui.getHandle('windowMain'), 'SelectionType');
      hListbox = this.gui.getHandle('listCalibrate');
     
      
      if strcmp(type, 'open')
        
        data = this.currentNode.handle.UserData.calibrate;
        
        cla(this.gui.getHandle('axesCalibrate'),'reset')
        axes(this.gui.getHandle('axesCalibrate'));
        hold all
        rotate3d on
        
        index = get(hListbox, 'Value');
        list = get(hListbox, 'String');
        
        for i = 1:numel(index)
          
          string = list{index(i)};
          C = textscan(string, '%s');
          index2 = cell2mat(data(:,4)) == str2double(C{1}{3});        
          
          track = cell2mat(data{index2, 1});
          
          if size(track,2) > 9
            plot3(track(:,7) .* 133, track(:,8) .* 133,track(:,11));

            view(3)
            %daspect([1,1,1])
          else
            waitfor(msgbox(['Track ' num2str(data{index2,4}) ' failed to calibrate'],...
                       'ERROR','error'));
          end
          
        end
        
        hold off
    
        axis equal
        axis vis3d
        box on
        
        
      end

end

% 
%     rotate3d on
%     
%     
%     for i = 1:n
%       
%       if ~selectedNodes(i).isRoot
%         data = selectedNodes(i).handle.UserData;
%         tracks = cell2mat(data{1,1});
%         
%         if size(tracks,2) > 8
%           plot3(tracks(:,7) .* 133, tracks(:,8) .* 133,tracks(:,11));
%           
%           view(3)
%           %daspect([1,1,1])
%         else
%           waitfor(msgbox(['Track ' num2str(data{1,4}) ' failed to calibrate'],...
%                      'ERROR','error'));
%         end
%       end
%       
%     end
%     
%     hold off
%     
%     axis equal
%     axis vis3d
%     box on

