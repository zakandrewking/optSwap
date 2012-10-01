function loopOptSwapD

    % setup Cleaner
    % cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'glyclate KO after swap';
    
    sets = [
        3,0;
        5,0;
           ];
    startSwaps = {
{'G6PDH2r';'LCARR';'TRSARr';}
{'GLUDy';'LCARR';'ME1';'PDH';}
                 };
        
    for i=1:size(sets,1)
        status = sprintf('run: %d, kos: %d swaps: %d', i, sets(i,1), sets(i,2)); 
        opt.knockoutNum = sets(i,1);
        opt.swapNum = sets(i,2);
        opt.targetRxns = {'EX_glyclt(e)'};
        opt.experiment = run;
        opt.logFile = 'database.csv';
        opt.startWithSwaps = startSwaps{i};
        runOptSwapD(opt);
    end
    status = 'finished';
end