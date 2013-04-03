function [newProb,solverParams] = setParams(prob, cobraSolverFlag, maxTime)

% OPTIONAL
% cobraSolverFlag - setup for solveCobraMILP
% maxTime - time limit in minutes (default 24 hours)
    if (nargin < 2 || isempty(cobraSolverFlag)), cobraSolverFlag = false; end
    if (nargin < 3 || isempty(maxTime)), maxTime = 24*60; end

    intTol = 10e-9;
    relMipGapTol = 1e-6;
    timeLimit = maxTime * 60; %seconds
    logFile = 'log.txt';
    printLevel = 10;
    feasTol = 1e-8;
    optTol = 1e-8;
    NUMERICALEMPHASIS = 1;
    absMipGapTol = 1e-8;

    if cobraSolverFlag
        newProb = prob;
        solverParams.relMipGapTol = relMipGapTol;
        solverParams.timeLimit = timeLimit;
        solverParams.logFile = logFile;
        solverParams.printLevel = printLevel;
        % if USE_MIP_FOCUS_1, solverParams.MIPFocus = 1;
        % else solverParams.MIPFocus = 3;
        % end
        % % solverParams.printLevel = printLevel;
        % solverParams.intTol = intTol;
        % solverParams.feasTol = feasTol;
        % solverParams.optTol = optTol;
        % solverParams.absMipGapTol = absMipGapTol;
        %gurobi default iteration limit is infinity
        
    else
        newProb=prob;
        newProb.Solver.Alg = 2; % Depth First, then Breadth (Default Depth
                                % First)
                                % 2: Depth first. When integer solution found, switch to Breadth.
        newProb.optParam.MaxIter = 100000; %100000
                                           % Must increase iterations from default 500
        newProb.optParam.IterPrint = 0;
        newProb.PriLev = printLevel;
        newProb.MIP.cpxControl.EPINT=intTol;
        newProb.MIP.cpxControl.EPOPT=optTol;
        newProb.MIP.cpxControl.EPRHS=feasTol;
        newProb.MIP.cpxControl.EPGAP=relMipGapTol;
        newProb.MIP.cpxControl.EPAGAP=absMipGapTol;
        newProb.MIP.cpxControl.NUMERICALEMPHASIS = ...
            NUMERICALEMPHASIS;
        newProb.MIP.cpxControl.SCAIND=1;
        newProb.MIP.cpxControl.TILIM=timeLimit;
    end
end