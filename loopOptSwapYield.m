function loopOptSwapYield
        
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'optSwapYield';
    logFile = sprintf('optSwapYield-Native-xylose_%s.tsv', ...
                      datestr(now, 'yy-mm-dd_HH_MM_SS'));
    global fileId
    fileId = fopen(logFile, 'a');
    fprintf(fileId, ['target\taerobic\tsubstrate\tnum swaps\tthko\' ...
                     'tf_k\tmax yield\tswaps\ttime (s)\n']);
    fclose(fileId);
    
    substrates = {'EX_xyl-D(e)'};
    aer = {'anaerobic','aerobic'};
    swaps = [0, 1, 2];
    for i=1:length(substrates)
        for j=1:length(aer)
            for k=1:length(swaps)
                opt.thko = 'nothko';
                opt.substrate = substrates{i};
                opt.aerobicString = aer{j};
                opt.swapNum = swaps(k);
                opt.logFile = logFile;
                runOptSwapYield(opt);
            end
        end
    end
end