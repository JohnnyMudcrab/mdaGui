% TODO: toolbar für figure
% TODO: consistenz vond aten sicher stellen (infobar)


% create mda object
hMda = mda();

% init gui
hMda.init();

% save object in appdata
setappdata(0,'hMda',hMda);

% change closing function in order to purge appdata correctly
set(gcf,'CloseRequestFcn',@hMda.closingFcn);

% % create toolbar
gui.initToolbar;

% set gui state to ready
hMda.gui.changeStatus('statusMain', 'Ready...');

% clear object from workspace
clear hMda