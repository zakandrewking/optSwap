function loopOptSwapDed

% setup Cleaner
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = '5 knocks on 3-swap acetate design';

    targetRxns = 'EX_ac(e)';
    sets = [5,0;];
    startWithSwaps = {'GLUDy';'MDH';'PGCD'};
    % % now convolute
    % sets = [];
    % for i=1:length(targetRxns), sets = [sets; setsOne]; end
    % fprintf('size of sets: %dx%d\n', size(sets,1), size(sets,2));
    for i=1:size(sets,1)
        status = sprintf('run: %d, kos: %d swaps: %d', i, sets(i,1), sets(i,2));
        opt.knockoutNum = sets(i,1);
        opt.swapNum = sets(i,2);
        opt.targetRxns = targetRxns;
        opt.experiment = run;
        opt.logFile = 'database-2.csv';
        opt.startWithSwaps = startWithSwaps;
        % opt.startWithKnocks = startKnocks;
        opt.swapAllDhs = false;
        runOptSwapD(opt);
    end
    status = 'finished';
end