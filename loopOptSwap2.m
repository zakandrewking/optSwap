function loopOptSwap2

% setup Cleaner
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';    
    run = 'best of 4--13 products--anaerobic aerobic xylose--dh kos ok';
    
    interventionNum = 4;
    aer = {'aerobic', 'anaerobic'};;
    substrates = {'EX_xyl-D(e)','EX_xyl-D(e)'};
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
        opt.useCobraSolver = true; 
        opt.allowDehydrogenaseKnockout = true;
        opt.logFile = 'database-2.csv';
        runOptSwap(opt);
    end
    status = 'finished';
end