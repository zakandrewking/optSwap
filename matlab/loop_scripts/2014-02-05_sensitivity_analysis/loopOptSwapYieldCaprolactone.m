function loopOptSwapYieldCaprolactone
    
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'optSwapYield';

    logFile = 'optSwapYield-caprolactone-aerobic.tsv';

    global fileId
    fileId = fopen(logFile, 'a');

    fprintf(fileId, ['target\taerobic\tsubstrate\tnum swaps\tthko\' ...
                     'tf_k\tmax yield\tswaps\ttime (s)\tsur\tbiomass\n']);
    fclose(fileId);
    
    max_sur = 20;
    count = 10;
    aer = 'aerobic';
    thko = 'nothko';
    substrate = 'EX_glc(e)';
    targetRxns = {'CMHO'};

    model = setupModel('iJO', substrate, aer, thko);
    model = changeRxnBounds(model, substrate, -max_sur, 'l');
    soln = optimizeCbModel(model);
    max_gr = soln.f;

    swaps = [0, 1, 2];
    minBiomasses = linspace(0, max_gr, count);
    surs = linspace(max_sur/count, max_sur, count)
    for i=1:length(minBiomasses)
        for j=1:length(surs)
            for k=1:length(swaps)
                opt.thko = thko;
                opt.substrate = substrate;            
                opt.aerobicString = aer;
                opt.logFile = logFile;
                opt.targetRxns = targetRxns;
                opt.caprolactone = true;

                opt.minBiomass = minBiomasses(i);
                opt.sur = surs(j);
                opt.swapNum = swaps(k);

                runOptSwapYield(opt);
            end
        end
    end
end