function buttonCalibrate(this)

  this.changeStatus('statusMain', 'Busy...');
  
  this.disableGui();
  set(this.getHandle('menuStop'), 'enable', 'on')

  hMda = getappdata(0, 'hMda');
  
  try
    hMda.calibrate();
  catch exc
    msgbox(getReport(exc), 'ERROR', 'error')
  end
  
  hMda.update
  
  this.enableGui();
  
  this.changeStatus('statusMain', 'Ready...');
  
end

