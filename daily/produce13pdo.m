addpath_recurse('~/repos/optSwap',{'.svn','.git'});

model = setupModel('iJO-h','EX_glyc(e)','anaerobic','thko');
model = turnOnSubSystem(model, 1, {'1,3-Propanediol production'});
printRxnFormula(model,model.rxns(ismember(model.subSystems,' 1,3-Propanediol production')));
% EX_13ppd(e)