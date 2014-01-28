function testGurobi

    run = 'test-gurobi';

    interventionNum = 1;
    aer = {'anaerobic'};
    substrates = {'EX_glc(e)'};
    i = 1;
    opt.knockoutNum = -1;
    opt.swapNum = -1;
    opt.interventionNum = interventionNum;
    opt.targetRxns = {'Ex_for(e)';};
    opt.experiment = run;
    opt.aerobicString = aer{i};
    opt.substrate = substrates{i};
    opt.maxTime = 1*60; %min
    opt.useCobraSolver = true;
    opt.allowDehydrogenaseKnockout = true;
    opt.logFile = 'database-test.csv';
    runOptSwap(opt);
end
