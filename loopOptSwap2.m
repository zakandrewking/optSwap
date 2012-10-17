function loopOptSwap2

% setup Cleaner
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'swaps and kos, EX_15dap(e)';

    targetRxns = 'EX_15dap(e)';
    sets = [
            3,3;
            5,5;
            ];
    for i=1:size(sets,1)
        status = sprintf('run: %d, kos: %d swaps: %d', i, sets(i,1), sets(i,2));
        opt.knockoutNum = sets(i,1);
        opt.swapNum = sets(i,2);
        opt.targetRxns = targetRxns;
        opt.experiment = run;
        opt.logFile = 'database-2.csv';
        % opt.startWithSwaps = startWithSwaps;
        % opt.startWithKnocks = startWithKnocks;
        opt.swapAllDhs = false;
        runOptSwapD(opt);
    end
    status = 'finished';
end