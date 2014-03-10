function loopOptSwap3

% setup Cleaner
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';    
    run = 'ethanol-glucose-anaerobic-1234-LX_OptSwap_LE';

    interventionNum = [-1, -1, -1, -1,  1,  2,  3,  4];
    knockoutNum     = [ 0,  0,  0,  0, -1, -1, -1, -1];
    swapNum =         [ 1,  2,  3,  4, -1, -1, -1, -1];
    
    aer = {'anaerobic'};
    substrates = {'EX_glc(e)'};
    for i=1:8
        status = sprintf('run %d', i);
        opt.knockoutNum = knockoutNum(i);
        opt.swapNum = swapNum(i);
        opt.interventionNum = interventionNum(i);
        opt.targetRxns = {'EX_etoh(e)'};
        opt.experiment = run;
        opt.aerobicString = aer{1};
        opt.substrate = substrates{1};
        opt.solverParams.maxTime = 12*60*60; %sec
        opt.solverParams.intTol = 1e-09;
        opt.solverParams.EPRHS = 1e-07;
        opt.solverParams.THREADS = 10;
        opt.useCobraSolver = false; 
        opt.allowDehydrogenaseKnockout = true;
        opt.logFile = 'database-3.csv';
        runOptSwap(opt);
    end
    status = 'finished';
end