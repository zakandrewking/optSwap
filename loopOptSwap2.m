function loopOptSwap2

% setup Cleaner
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';    
    run = 'ethanol-glucose-anaerobic-test cobraSolverFlag';

    aer = {'anaerobic'};
    substrates = {'EX_glc(e)'};
    for i=1:2
        status = sprintf('run %d', i);
        opt.knockoutNum = -1;
        opt.swapNum = -1;
        opt.interventionNum = 1;
        opt.targetRxns = {'EX_etoh(e)'};
        opt.experiment = run;
        opt.aerobicString = aer{1};
        opt.substrate = substrates{1};
        opt.solverParams.maxTime = 12*60*60; %sec
        opt.solverParams.intTol = 1e-09;
        opt.solverParams.EPRHS = 1e-07;
        opt.solverParams.THREADS = 10;
        opt.useCobraSolver = false;
        opt.allowDehydrogenaseKnockout = false;
        opt.logFile = 'database-2.csv';
        if i==1
            opt.useCobraSolver = true;
        else
            opt.useCobraSolver = false;
        end
        runOptSwap(opt);
    end
    status = 'finished';
end