function loopOptSwapYieldLycopene
        
    run = 'optSwapYieldLycopene';
    logFile = sprintf('optSwapYield-lycopene_%s.tsv', ...
                      datestr(now, 'yy-mm-dd_HH_MM_SS'));
    global fileId
    fileId = fopen(logFile, 'a');
    fprintf(fileId, ['target\tsubsystem\taerobic\tsubstrate\tnum swaps\tthko\' ...
                     'tf_k\tmax yield\tswaps\ttime (s)\n']);
    fclose(fileId);

    substrates = {'EX_glc(e)', 'EX_xyl_D(e)', 'EX_glyc(e)'};
    aer = {'anaerobic','aerobic'};
    swaps = [0, 1, 2];
    for i=1:length(substrates)
        for j=1:length(aer)
            model = setupModel('iJO',substrates{i},aer{j},'nothko');
            model = makeLycopene(model)
            soln = optimizeCbModel(model);
            if (i==3 && j==1) %glycerol anaerobic
                opt.minBiomass = 0.1*soln.f;
            else
                opt.minBiomass = 0.1;
            end
            fprintf('minBiomass\t%.3f', opt.minBiomass);
            for k=1:length(swaps)                
                opt.thko = 'nothko';
                opt.substrate = substrates{i};
                opt.aerobicString = aer{j};
                opt.swapNum = swaps(k);
                opt.logFile = logFile;
                opt.lycopene = true;
                opt.targetRxns = {'EX_lyco(e)'};
                runOptSwapYield(opt);
            end
        end
    end
end