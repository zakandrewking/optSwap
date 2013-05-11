function [newProb,solverParams] = setParams(prob, cobraSolverFlag, nonDefaultParameters)

% OPTIONAL
% cobraSolverFlag - setup for solveCobraMILP
% maxTime - time limit in minutes (default 24 hours)

    if (nargin < 2 || isempty(cobraSolverFlag)), cobraSolverFlag = false; end

    intTol = 10e-12;
    relMipGapTol = 1e-6;
    logFile = 'log.txt';
    printLevel = 10;
    feasTol = 1e-8;
    optTol = 1e-8;
    NUMERICALEMPHASIS = 1;
    absMipGapTol = 1e-8;
    EPRHS = 1e-9;
    maxTime = 12*60*60; % seconds
    THREADS = 10;
    
    if isfield(nonDefaultParameters, 'maxTime')
        maxTime = nonDefaultParameters.maxTime;
        fprintf('max time: %d seconds\n', maxTime);
    end
    if isfield(nonDefaultParameters, 'intTol')
        intTol = nonDefaultParameters.intTol;
        fprintf('intTol: %.1e\n', intTol);
    end
    if isfield(nonDefaultParameters, 'EPRHS')
        EPRHS = nonDefaultParameters.EPRHS;
        fprintf('EPRHS: %.1e\n', EPRHS);
    end
    if isfield(nonDefaultParameters, 'THREADS')
        THREADS = nonDefaultParameters.THREADS;
        fprintf('THREADS: %d\n', THREADS);
    end
    if cobraSolverFlag
        newProb = prob;
        solverParams.relMipGapTol = relMipGapTol;
        solverParams.timeLimit = maxTime;
        solverParams.logFile = logFile;
        solverParams.printLevel = printLevel;
        solverParams.intTol = intTol;
        solverParams.feasTol = feasTol;
        solverParams.optTol = optTol;
        solverParams.absMipGapTol = absMipGapTol;
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
        newProb.MIP.cpxControl.TILIM=maxTime;
        newProb.MIP.cpxControl.EPRHS = EPRHS;
        newProb.MIP.cpxControl.THREADS = THREADS;
        newProb.CPLEX.LogFile = 'cplex-log.txt';
    end
end