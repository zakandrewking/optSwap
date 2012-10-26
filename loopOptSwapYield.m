function loopOptSwapYield
        
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'optSwapYield--all one and two swaps';
    
    substrates = {'EX_glc(e)', 'EX_xyl-D(e)'};
    aer = {'anaerobic','aerobic'};
    swaps = [1,2];
    for i=1:length(substrates)
        for j=1:length(aer)
            for k=1:length(swaps)
                opt.substrate = substrates{i};
                opt.aerobicString = aer{j};
                opt.swapNum = swaps(k);
                opt.logFile = 'optSwapYield-database.csv';
                runOptSwapYield(opt);
            end
        end
    end
end