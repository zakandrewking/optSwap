function loopOptSwapWithSwaps

%setup Cleaner 
cleaner = onCleanup(@() cleanup);
global run status
status = 'starting';
run = '3 knocks after 4 swaps';
    
swaps = {
    {'EAR60x';'GLUDy';'ME1';'TRSARr';}
    {'3OAR40';'G6PDH2r';'GLUDy';'TRSARr';}
    {'G6PDH2r';'MDH';'MTHFD';'PDH';}
    {'GLUDy';'MDH';'PDH';'PGCD';}
        };  
targetRxns = {
    'EX_etoh(e)'
    'EX_for(e)'
    'EX_succ(e)'
    'EX_ac(e)'
             };
koNum = 3;

for i = 1:length(targetRxns)
    status = sprintf('%d: %s', i, targetRxns{i});
    opt.startWithSwaps = swaps{i};
    opt.targetRxns = targetRxns{i};
    opt.knockoutNum = koNum;
    opt.experiment = run;
    swapStr = '';
    for j = 1:length(swaps{i}), swapStr = [swapStr ' ' swaps{i}{j}]; end
    opt.notes = swapStr;
    runOptSwapD(opt);
end

status = 'finished';
end