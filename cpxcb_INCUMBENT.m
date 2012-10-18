%function ret = cpxcb_INCUMBENT(x,f,Prob,cbxCBInfo)
%
% CPLEX MIP Incumbent callback
%
% Called from TOMLAB /CPLEX during mixed integer optimization when a new integer
% solution has been found but before this solution has replaced the current best known integer solution.
%
% This file can be used to perform any desired analysis of the new integer
% solution and return a status flag to the solver deciding whether to stop
% or continue the optimization, and also whether to accept or discard the newly
% found solution.
%
% This callback is enabled by setting callback(14)=1 in the call to
% cplex.m, or Prob.MIP.callback(14)=1 if using tomRun('cplex',...)
%
% cpxcb_INCUMBENT is called by the solver with three arguments:
%
%  x    - the new integer solution
%  f    - the objective value at x
%  Prob - the Tomlab problem structure
%
% cpxcb_INCUMBENT should return one of the following scalar values:
%
%   0    Continue optimization and accept new integer solution
%   1    Continue optimization but discard new integer solution
%   2    Stop optimization and accept new integer solution
%   3    Stop optimization adn discard new integer solution
%
% Any other return value will be interpreted as 0.
%
% If modifying this file, it is recommended to make a copy of it which
% is placed before the original file in the MATLAB path.
%

% Anders Goran, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2007 by Tomlab Optimization Inc., $Release: 10.1.0$
% Written Jun 1, 2007.  Last modified Jun 1, 2007.

function ret = cpxcb_INCUMBENT(x,f,Prob)

% ADD USER CODE HERE.

% Accepted return values are:
%
%   0    Continue optimization and accept new integer solution
%   1    Continue optimization but discard new integer solution
%   2    Stop optimization and accept new integer solution
%   3    Stop optimization and discard new integer solution
%
% Any other return value will be interpreted as 0.

global optSwapCallbackOptions
global solID

file = optSwapCallbackOptions.intermediateSolutionsFile;
global fileId
fileId = fopen(intermediateSolutionsFile, 'a');

% Initialize
if isempty(solID)
    solID = 0;
    objective = [];
    growth = [];
    koList = {};
    swapList = {};
    myPrint('%s\n', intermediateSolutionsFile);
    myPrint('time/date,f,kos,swaps', []);
end

solID = solID + 1;

% Get the reactions
model = optSwapCallbackOptions.model;
yInd = optSwapCallbackOptions.yInd;
qInd = optSwapCallbackOptions.qInd;
intVars = optSwapCallbackOptions.intVars;

y = x(intVars(1:length(yInd)));
koRxns = model.rxns(yInd(y==0));
q = x(intVars(length(yInd)+1:length(yInd)+length(qInd)));
swapRxns = model.rxns(qInd(q==0));

myPrint('\n', [])
myPrint('%s,', datestr(now));
myPrint('%.4f,', f);
printReactions(koRxns);
printReactions(swapList);
ret = 0;
fclose(fileId);
end

function printReactions(reactions)
for j=1:length(reactions)
    myPrint('''%s'';', reactions{j});
end
myPrint(',', []);
end

function myPrint(string, val)
global fileId
display(sprintf(string, val));
fprintf(fileId, string, val);
end
