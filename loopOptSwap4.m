function loopOptSwap4

% setup Cleaner
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';    
    run = 'best of 1--all targets with new dhlist - glc';
    
    interventionNum = 1;
    aer = {'anaerobic','aerobic','anaerobic','aerobic'};
    substrates = {'EX_glc(e)','EX_glc(e)'};
    for i=1:2
        status = sprintf('run %d: %d intervention(s)', i, interventionNum);
        opt.knockoutNum = -1;
        opt.swapNum = -1;
        opt.interventionNum = interventionNum;
        opt.targetRxns = {'EX_etoh(e)';
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
                          'EX_glu-L(e)';
                         };
        opt.experiment = run;
        opt.aerobicString = aer{i};
        opt.substrate = substrates{i};
        opt.maxTime = 12*60; %min
        opt.useCobraSolver = false;
        opt.logFile = 'database-4.csv';
        runOptSwap(opt);
    end
    status = 'finished';
end