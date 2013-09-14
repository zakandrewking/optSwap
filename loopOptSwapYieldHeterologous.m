function loopOptSwapYieldHeterologous
        
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'optSwapYield';
    logFile = sprintf('optSwapYield-13pdo_%s.tsv', ...
                      datestr(now, 'yy-mm-dd_HH_MM_SS'));
    global fileId
    fileId = fopen(logFile, 'a');
    fprintf(fileId, ['target\tsubsystem\taerobic\tsubstrate\tnum swaps\tthko\' ...
                     'tf_k\tmax yield\tswaps\ttime (s)\n']);
    fclose(fileId);
    
    substrates = {'EX_glc(e)', 'EX_xyl_D(e)', 'EX_glyc(e)'};
    aer = {'anaerobic','aerobic'};
    swaps = [0, 0, 1, 2];
    targetRxns = {'EX_13ppd(e)'};
    for i=1:length(substrates)
        for j=1:length(aer)
            for k=1:length(swaps)
                if k==1
                    opt.thko = 'nothko';
                else
                    opt.thko = 'thko';
                end
                opt.substrate = substrates{i};
                opt.aerobicString = aer{j};
                opt.swapNum = swaps(k);
                opt.targetRxns = targetRxns;
                opt.logFile = logFile;
                runOptSwapYieldHeterologous(opt);
            end
        end
    end
end