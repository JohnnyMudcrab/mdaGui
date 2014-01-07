function buttonLocate(this)
  
this.changeStatus('statusMain', 'Busy...');
  %this.disableGui();

  hMda = getappdata(0, 'hMda');
  hMda.locate(true);
  
  %this.enableGui();
  this.changeStatus('statusMain', 'Ready...');
  
end

