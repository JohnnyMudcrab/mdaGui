function buttonTrack(this)

  this.changeStatus('statusMain', 'Busy...');
  this.disableGui();
  set(this.getHandle('menuStop'), 'enable', 'on')

  hMda = getappdata(0, 'hMda');
  hMda.track(true);
  
  this.enableGui();
  this.changeStatus('statusMain', 'Ready...');
  
end
