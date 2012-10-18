function loopOptSwap1

    % setup Cleaner
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'TEST intermediate solns, 3x3 on ethanol';

    sets = [2,2];

    for i=1:size(sets,1)
        status = sprintf('run: %d, kos: %d swaps: %d', i, sets(i,1), sets(i,2));
        opt.knockoutNum = sets(i,1);
        opt.swapNum = sets(i,2);
        opt.targetRxns = {'EX_etoh(e)'};
        opt.experiment = run;
        opt.logFile = 'database-1.csv';
        opt.aerobicString = 'anaerobic';
        opt.substrate = 'EX_glc(e)';
        opt.maxTime = 168*60; %min
        opt.swapAllDhs = false;
        opt.printIntermediateSolutions = true;
        runOptSwap(opt);
    end
    status = 'finished';
end