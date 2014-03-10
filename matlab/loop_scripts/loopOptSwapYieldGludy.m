function loopOptSwapYieldGludy
        
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'optSwapYield';
    logFile = sprintf('optSwapYield_GLUDy-out_%s.tsv', ...
                      datestr(now, 'yy-mm-dd_HH_MM_SS'));
    substrates = {'EX_glc(e)', 'EX_xyl-D(e)', 'EX_glyc(e)'};
    aer = {'anaerobic','aerobic'};
    swaps = [1, 2, 3];
    for i=1:length(substrates)
        for j=1:length(aer)
            for k=1:length(swaps)
                opt.substrate = substrates{i};
                opt.aerobicString = aer{j};
                opt.swapNum = swaps(k);
                opt.logFile = logFile;
                opt.dhRxns = dhRxnList('yield-gludy-out')
                runOptSwapYield(opt);
            end
        end
    end
end