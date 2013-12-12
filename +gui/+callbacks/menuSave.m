function menuSave(this)

  this.changeStatus('statusMain', 'Busy...');

  % get tree
  hTree = this.getHandle('treeMain');
  
  % save tree
  hTree.save(); 
  
  this.changeStatus('statusMain', 'Ready...');

end

