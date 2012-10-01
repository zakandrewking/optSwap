function loopOptSwapDed

    % setup Cleaner
    % cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'glyclate production';
    
    sets = [
        0,3;
        0,5;
           ];
%     startSwaps = {
% {'G6PDH2r';'GLUDy';'MTHFD';}
% {'3OAR60';'3OAR80';'AKGDH';'GLUDy';}
%         };
        
    for i=1:size(sets,1)
        status = sprintf('run: %d, kos: %d swaps: %d', i, sets(i,1), sets(i,2)); 
        opt.knockoutNum = sets(i,1);
        opt.swapNum = sets(i,2);
        opt.targetRxns = {'EX_glyclt(e)'};
        opt.experiment = run;
        opt.logFile = 'database-2.csv';
        % opt.startWithSwaps = startSwaps{i};
        runOptSwapD(opt);
    end
    status = 'finished';
end