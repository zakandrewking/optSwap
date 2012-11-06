function loopOptSwap1

% setup Cleaner
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'TEST RobustKnock flag';

    % sets = [2,0;
    %         2,0;
    %         4,0;
    %         4,0];
    % aer = {'anaerobic',
    %        'aerobic',
    %        'anaerobic',
    %        'aerobic'};
    sets = [4,0];
    for i=1:1
        status = sprintf('run: %d, kos: %d swaps: %d', i, sets(i,1), sets(i,2));
        opt.knockoutNum = sets(i,1);
        opt.swapNum = sets(i,2);
        opt.targetRxns = 'EX_etoh(e)';
                          % 'EX_for(e)';
                          % 'EX_succ(e)';
                          % 'EX_ac(e)';
                          % 'EX_lac-D(e)';
                          % 'EX_akg(e)';
                          % 'EX_ala-L(e)';
                          % 'EX_glyc(e)';
                          % 'EX_ser-L(e)';
                          % 'EX_pyr(e)';
                          % 'EX_fum(e)';
                          % 'EX_mal-L(e)';
                          % 'EX_glu-L(e)';}; 
        opt.experiment = run;
        opt.logFile = 'database-1.csv';
        opt.aerobicString = 'anaerobic';
        opt.substrate = 'EX_glc(e)';
        opt.maxTime = 12*60; %min
        opt.useCobraSolver = true;
        opt.canKnockDHs = true; 
        opt.knockType = 2
        runOptSwap(opt);
    end
    status = 'finished';
end