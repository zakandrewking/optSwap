function loopOptSwap1

% setup Cleaner
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';    
    
    runs = {
        {4},{'EX_for(e)'},{'anaerobic'},{'EX_glc(e)'};
        {4},{'EX_akg(e)'},{'anaerobic'},{'EX_glc(e)'};
        {4},{'EX_for(e)'},{'anaerobic'},{'EX_xyl-D(e)'};
        {4},{'EX_succ(e)'},{'anaerobic'},{'EX_xyl-D(e)'};
        {4},{'EX_akg(e)'},{'anaerobic'},{'EX_xyl-D(e)'};
        {4},{'EX_etoh(e)'},{'aerobic'},{'EX_glc(e)'};
        {4},{'EX_for(e)'},{'aerobic'},{'EX_glc(e)'};
        {4},{'EX_succ(e)'},{'aerobic'},{'EX_glc(e)'};
        {4},{'EX_ac(e)'},{'aerobic'},{'EX_glc(e)'};
        {4},{'EX_lac-D(e)'},{'aerobic'},{'EX_glc(e)'};
        {4},{'EX_etoh(e)'},{'anaerobic'},{'EX_xyl-D(e)'};
        {4},{'EX_for(e)'},{'aerobic'},{'EX_glc(e)'}
           };
    run = 'rerunning best of 4s'
   
    for i=1:4
        status = sprintf('run %d', i);
        opt.knockoutNum = -1;
        opt.swapNum = -1;
        opt.interventionNum = runs{i,1}{1};
        opt.targetRxns = runs{i,2};
        opt.experiment = run;
        opt.aerobicString = runs{i,3}{1};
        opt.substrate = runs{i,4}{1};
        opt.solverParams.maxTime = 12*60*60; %sec
        opt.solverParams.intTol = 1e-09; %sec
        opt.solverParams.EPRHS = 1e-09; 
        opt.solverParams.PROBE = 3;
        opt.solverParams.THREADS = 20;
        opt.useCobraSolver = false; 
        opt.allowDehydrogenaseKnockout = true;
        opt.logFile = 'database-1.csv';
        runOptSwap(opt);
    end
    status = 'finished';
end