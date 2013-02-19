function loopOptSwapYield
        
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'optSwapYield';
    logFile = sprintf('optSwapYield_no_GAPD_%s.txt', ...
                      datestr(now, 'yy-mm-dd_HH_MM_SS'));
    substrates = {'EX_glc(e)', 'EX_xyl-D(e)'};
    aer = {'anaerobic','aerobic'};
    swaps = [1, 2];
    for i=1:length(substrates)
        for j=1:length(aer)
            for k=1:length(swaps)
                opt.substrate = substrates{i};
                opt.aerobicString = aer{j};
                opt.swapNum = swaps(k);
                opt.logFile = logFile;
                opt.dhRxns = dhRxnList(20);
                runOptSwapYield(opt);
            end
        end
    end
end