function loopOptSwap1

% setup Cleaner
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = '1,5dap one and two swaps';

    sets = [0,3;
            0,4];
    for i=1:size(sets,1)
        status = sprintf('run: %d, kos: %d swaps: %d', i, sets(i,1), sets(i,2));
        opt.knockoutNum = sets(i,1);
        opt.swapNum = sets(i,2);
        opt.targetRxns = {'EX_15dap(e)'};
        opt.experiment = run;
        opt.logFile = 'database-1.csv';
        opt.aerobicString = 'anaerobic';
        opt.substrate = 'EX_glc(e)';
        opt.maxTime = 12*60; %min
        opt.useCobraSolver = true;
        runOptSwap(opt);
    end

    status = 'finished';
end