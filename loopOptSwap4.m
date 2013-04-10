function loopOptSwap4

% setup Cleaner
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';    
    run = ['RK--for formate anaerobic  glucose--dh kos ok--intTol1e-12 EPRHS1e-8'];
    
    aer = {'anaerobic', 'aerobic'};
    substrates = {'EX_glc(e)', 'EX_glc(e)'};
    for i=1:1
        status = sprintf('run %d: 3 knockouts(s)', i);
        opt.knockoutNum = 3;
        opt.swapNum = 0;
        opt.interventionNum = -1;
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
        opt.solverParams.maxTime = 12*60*60; %sec
        opt.useCobraSolver = false;
        opt.logFile = 'database-4.csv';
        runOptSwap(opt);
    end
    status = 'finished';
end