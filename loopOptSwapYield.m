function loopOptSwapYield
        
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'optSwapYield--all two swaps--no MDH';
    logFile = sprintf('optSwapYield_%s.txt', ...
                      datestr(now, 'yy-mm-dd_HH_MM_SS'));
    substrates = {'EX_glc(e)', 'EX_xyl-D(e)'};
    aer = {'anaerobic','aerobic'};
    swaps = [1, 2];
    for i=1:length(substrates)
        for j=1:length(aer)
            for k=1:length(swaps)
                opt.substrate = substrates{i};
                opt.aerobicString = aer{j};
                opt.swapNum = s waps(k);
                opt.logFile = logFile;
                opt.dhRxns = {'3OAR100';
                '3OAR40';
                '3OAR60';
                '3OAR80';
                'AGPR';
                'AKGDH';
                'ASAD';
                'DHDPRy';
                'EAR40x';
                'EAR60x';
                'EAR80x';
                'G6PDH2r';
                'GAPD';
                % 'GLUDy' removed because of loop result
                'GND';
                'HPYRRx';
                'HSDy';
                'ICDHyr';
                'IPMD';
                'KARA1';
                'KARA2';
                'MDH'; 
                'MTHFD';
                'NADH16pp';
                'PDH';
                'PGCD';
                'SHK3Dr';
                % 'TRSARr'; removed because of loop result
                'ME1';
                'LDH_D';
                'LCARR';};
                
                runOptSwapYield(opt);
            end
        end
    end
end