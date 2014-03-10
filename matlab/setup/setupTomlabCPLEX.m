addpath /Library/tomlab/tomsym
addpath /Library/tomlab/tomsym/examples
addpath /Library/tomlab/tomsym/funcs
addpath /Library/tomlab/tomsym/userfun_example
addpath /Library/tomlab/gurobi/examples
addpath /Library/tomlab/gurobi
addpath /Library/tomlab/cplex/network
addpath /Library/tomlab/cplex/examples
addpath /Library/tomlab/cplex
addpath /Library/tomlab
addpath /Library/tomlab/mad
addpath /Library/tomlab/modellib
addpath /Library/tomlab/modellib/airtransport
addpath /Library/tomlab/modellib/economics
addpath /Library/tomlab/modellib/groundtransport
addpath /Library/tomlab/modellib/loadingandcutting
addpath /Library/tomlab/modellib/miningandprocess
addpath /Library/tomlab/modellib/planning
addpath /Library/tomlab/modellib/publicservices
addpath /Library/tomlab/modellib/scheduling
addpath /Library/tomlab/modellib/telecommunication
addpath /Library/tomlab/modellib/timetabling
addpath /Library/tomlab/base
addpath /Library/tomlab/lib
addpath /Library/tomlab/cgo
addpath /Library/tomlab/mex
addpath /Library/tomlab/testprob
addpath /Library/tomlab/examples
addpath /Library/tomlab/quickguide
addpath /Library/tomlab/usersguide
addpath /Library/tomlab/ampl
addpath /Library/tomlab/common
addpath /Library/tomlab/finance
addpath /Library/tomlab/optim
solverOK = false;
if ~solverOK && exist('tomRun.m','file')
    solverOK = changeCobraSolver('tomlab_cplex','LP');
    solver = 'tomlab_cplex';
end
if solverOK
    display(['Setting up solver ' solver]);
    changeCobraSolver(solver,'MILP');
    changeCobraSolver(solver,'QP');
    changeCobraSolver(solver,'MIQP');
else
    error('no solver available');
end