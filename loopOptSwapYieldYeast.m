function loopOptSwapYieldYeast
        
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'optSwapYieldYeast';
    logFile = sprintf('optSwapYieldYeast_%s.tsv', ...
                      datestr(now, 'yy-mm-dd_HH_MM_SS'));
    global fileId
    fileId = fopen(logFile, 'a');
    fprintf(fileId, ['target\taerobic\tsubstrate\tnum swaps\tthko\' ...
                     'tf_k\tmax yield\tswaps\ttime (s)\n']);
    fclose(fileId);
    
    % glc, D-xyl, gylc, L-arab
    substrates = {'r_1714', 'r_1718', 'r_1808', 'r_1878'};
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
                opt.modelname = 'yeast6';
                opt.dhRxns = yeastDhPool();
                runOptSwapYield(opt);
            end
        end
    end
end

function pool = yeastDhPool()
    pool = {}
end