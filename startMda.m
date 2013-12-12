% TODO: zoom in figure für all frames gleich
% TODO: was wenn current image entfernt wird?
% TODO: consistenz von daten sicher stellen (infobar)


% create mda object
hMda = mda();

% init gui
hMda.init();

% save object in appdata
setappdata(0,'hMda',hMda);

% change closing function in order to purge appdata correctly
set(gcf,'CloseRequestFcn',@hMda.closingFcn);

% create toolbar
gui.initToolbar(hMda.gui);

% set gui state to ready
hMda.gui.changeStatus('statusMain', 'Ready...');

% clear object from workspace
clear hMda