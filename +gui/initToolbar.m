function initToolbar(this)

  % enable figure toolbar
  set(gcf, 'toolbar', 'figure' )

  hToolbar = findall(gcf,'tag','FigureToolBar');

  load('+gui/icon.mat')
  load('+gui/icon2.mat')
  uipushtool(hToolbar,'CData',icon, 'Enable', 'on', 'ClickedCallback', @(src,event)gui.callbacks.toolbarApply(this))
  uipushtool(hToolbar,'CData',icon2, 'Enable', 'on', 'ClickedCallback', @(src,event)gui.callbacks.toolbarStop(this))
  
  for i = 1:25
    uipushtool(hToolbar,'CData',[], 'Enable', 'off')
  end

  % delete unnecessary toolbar icons and seperators
  set(findall(hToolbar,'ToolTipString','New Figure'), 'Visible','Off');
  set(findall(hToolbar,'ToolTipString','Open File'), 'Visible','Off');
  set(findall(hToolbar,'ToolTipString','Save Figure'), 'Visible','Off');
  set(findall(hToolbar,'ToolTipString','Print Figure'), 'Visible','Off');
  set(findall(hToolbar,'ToolTipString','Insert Legend'), 'Visible','Off');
  set(findall(hToolbar,'ToolTipString','Insert Colorbar'), 'Visible','Off');
  set(findall(hToolbar,'ToolTipString','Insert Colorbar'), 'Separator','Off');
  set(findall(hToolbar,'ToolTipString','Link Plot'), 'Visible','Off');
  set(findall(hToolbar,'ToolTipString','Link Plot'), 'Separator','Off');
  set(findall(hToolbar,'ToolTipString','Hide Plot Tools'), 'Visible','Off');
  set(findall(hToolbar,'ToolTipString','Show Plot Tools and Dock Figure'), 'Visible','Off');
  set(findall(hToolbar,'ToolTipString','Show Plot Tools and Dock Figure'), 'Separator','Off');

  % change order of pushtools
  oldOrder = allchild(hToolbar);
  newOrder = [oldOrder(28:end); oldOrder(1:27)];
  set(hToolbar,'Children',newOrder);

end