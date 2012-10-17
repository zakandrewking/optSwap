function loopOptSwap3

% setup Cleaner
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = '1x1 on all 13 important products-Aerobic';

    targetRxns = {'EX_etoh(e)';
                  'EX_for(e)';
                  'EX_succ(e)';
                  'EX_ac(e)';
                  'EX_lac-D(e)';
                  'EX_akg(e)';
                  'EX_ala-L(e)';
                  'EX_glyc(e)';
                  'EX_ser-L(e)';
                  'EX_pyr(e)';
                  'EX_fum(e)';
                  'EX_mal-L(e)';
                  'EX_glu-L(e)';};
    sets = [1,1];
    for i=1:size(sets,1)
        status = sprintf('run: %d, kos: %d swaps: %d', i, sets(i,1), sets(i,2));
        opt.knockoutNum = sets(i,1);
        opt.swapNum = sets(i,2);
        opt.targetRxns = targetRxns;
        opt.experiment = run;
        opt.logFile = 'database-3.csv';
        opt.swapAllDhs = false;
        opt.aerobicString = 'aerobic';
        runOptSwapD(opt);
    end
    status = 'finished';
end