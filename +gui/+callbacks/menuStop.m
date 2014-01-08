function menuStop(this)
  
  this.changeStatus('statusMain', 'Busy...Preparing to Stop, Please Wait!');
  
  drawnow

  hMda = getappdata(0, 'hMda');
  
  hMda.stop = true;
  
end

