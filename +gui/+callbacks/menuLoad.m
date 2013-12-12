function menuLoad(this)

  this.changeStatus('statusMain', 'Busy...');
  
  % get tree
  hTree = this.getHandle('treeMain');
  
  % save tree
  hTree.load();
  
  this.changeStatus('statusMain', 'Ready...');
  
end

