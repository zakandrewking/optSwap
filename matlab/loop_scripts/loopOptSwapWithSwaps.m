function loopOptSwapWithSwaps

%setup Cleaner 
cleaner = onCleanup(@() cleanup);
global run status
status = 'starting';
run = '5 knocks after 5 swaps';
    
swaps = {
    {'3OAR40';'3OAR80';'GLUDy';'LCARR';'ME1';}
    {'3OAR40';'EAR60x';'G6PDH2r';'GLUDy';'TRSARr';}
    {'G6PDH2r';'LCARR';'MDH';'MTHFD';'PDH';}
    {'GLUDy';'GND';'LCARR';'MDH';'PGCD';}
        };  
targetRxns = {
    'EX_etoh(e)'
    'EX_for(e)'
    'EX_succ(e)'
    'EX_ac(e)'
             };
koNum = 5;

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