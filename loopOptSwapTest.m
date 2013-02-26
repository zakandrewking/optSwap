function loopOptSwapTest

% setup Cleaner
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';    
    run = 'test-master-3-kos-formate-759d8';
    
    interventionNum = 3;
    aer = {'aerobic'};
    substrates = {'EX_xyl-D(e)','EX_glc(e)'};
    for i=1:1
        status = sprintf('run %d: %d intervention(s)', i, interventionNum);
        opt.knockoutNum = -1;
        opt.swapNum = -1;
        opt.interventionNum = interventionNum;
        opt.targetRxns = {% 'EX_etoh(e)';
                          'EX_for(e)';
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
                          % 'EX_glu-L(e)';
                         };
        opt.experiment = run;
        opt.aerobicString = aer{i};
        opt.substrate = substrates{i};
        opt.maxTime = 12*60; %min
        opt.useCobraSolver = false; 
        opt.allowDehydrogenaseKnockout = true;
        opt.logFile = 'database-0.csv';
        runOptSwap(opt);
    end
    status = 'finished';
end
