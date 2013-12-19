function buttonLocateApplyAll(this)

  this.changeStatus('statusMain', 'Busy...');
  this.disableGui();

  hMda = getappdata(0, 'hMda');
  hMda.locate();
  
  this.enableGui();
  this.changeStatus('statusMain', 'Ready...');

end

