function loopOptSwapDed

    % setup Cleaner
    % cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'glycolate production';
    
    sets = [
        0,3;
        0,5;
        3,0;
        5,0;
        3,3;
        5,5;]
        
    for i=1:size(sets,1)
        status = sprintf('run: %d, kos: %d swaps: %d', i, sets(i,1), sets(i,2));
        % for j=1:length(swaps)
        opt.knockoutNum = sets(i,1);
        opt.swapNum = sets(i,2);
        opt.targetRxns = {'EX_ptrc(e)'};
        opt.experiment = run;
        opt.logFile = 'database-2.csv';
        runOptSwapDed(opt);
        % end
    end
    status = 'finished'
end