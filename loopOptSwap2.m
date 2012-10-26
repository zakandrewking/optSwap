function loopOptSwap2

% setup Cleaner
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'TEST gurobi5 def. params--2x2 ethanol';

    sets = [2,2];
    global USE_MIP_FOCUS_1
    USE_MIP_FOCUS_1 = false;

    status = sprintf('run: %d, kos: %d swaps: %d', i, sets(1), sets(2));
    opt.knockoutNum = sets(1);
    opt.swapNum = sets(2);
    opt.targetRxns = {'EX_etoh(e)'};
    opt.experiment = run;
    opt.logFile = 'database-2.csv';
    opt.aerobicString = 'anaerobic';
    opt.substrate = 'EX_glc(e)';
    opt.maxTime = 12*60; %min
    opt.useCobraSolver = true;
    opt.notes = 'MIPFocus = 3';
    runOptSwap(opt);

    status = 'finished';
end