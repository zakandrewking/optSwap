function loopOptSwapMinRxns

    % setup Cleaner
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'Speed of KOs with 20 rxn set';

    sets = [
        1,0;
        2,0;
        3,0;
        5,0;
        10,0;
        20,0;
           ];

    for i=1:size(sets,1)
        status = sprintf('run: %d, kos: %d swaps: %d', i, sets(i,1), sets(i,2));
        opt.knockoutNum = sets(i,1);
        opt.swapNum = sets(i,2);
        opt.targetRxns = {'EX_for(e)';
                          'EX_ac(e)';};
        opt.experiment = run;
        opt.logFile = 'database.csv';
        reactionSet = koRxnList(20);
        runOptSwapRxnSet(opt, reactionSet);
    end
    status = 'finished';

end