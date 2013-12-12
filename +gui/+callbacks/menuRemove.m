function menuRemove(this)

  this.changeStatus('statusMain', 'Busy...');

  % get tree
  hTree = this.getHandle('treeMain');
  
  % save tree
  hTree.removeSelectedNodes(); 
  
  this.changeStatus('statusMain', 'Ready...');
  
end