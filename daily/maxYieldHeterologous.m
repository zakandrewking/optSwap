function maxYieldHeterologous()
    logFile = sprintf('maxYields-Heterologous_%s.tsv', ...
                       datestr(now, 'yy-mm-dd_HH_MM_SS'));
    global fileId
    fileId = fopen(logFile, 'a');
    fprintf(fileId, ['target\tsubsystem\taerobic\tsubstrate\tthko\tyield\n']);

    subs = {'EX_glc(e)','EX_xyl_D(e)'};
    aer = {'aerobic','anaerobic'};
    thko = {'thko','nothko'};
    targets = returnTargetRxns();
    subSystems = returnSubSystems();
    for i = 1:length(subs)
        for j = 1:length(aer)
            for k = 1:length(thko)
                [model, biomass] = setupModel('iJO-h',subs{i},aer{j},thko{k});
                model = changeRxnBounds(model,biomass,0.1,'l');
                for m = 1:length(targets)
                    modelT = model;
                    modelT = changeObjective(modelT,targets{m});
                    modelT = turnOnSubSystem(modelT, m, subSystems);
                    soln = optimizeCbModel(modelT);
                    myPrint('%s\t',targets{m});
                    myPrint('%s\t', subSystems{i});
                    myPrint('%s\t', aer{j});
                    myPrint('%s\t', subs{i});
                    myPrint('%s\t', thko{k});
                    myPrint('%.4f\t', soln.f);
                    myPrint('\n',[]);
                end
            end
        end
    end

    fclose(fileId);
end

function myPrint(string, val)
    global fileId
    fprintf(string, val);
    fprintf(fileId, string, val);
end

function rxns = returnTargetRxns()
    rxns = {'EX_13ppd(e)';
            'EX_14btd(e)';
            'EX_btd_RR(e)';
            'EX_btd_SS(e)';
            'EX_btd_meso(e)';
            'EX_bhb(e)';
            'EX_bhb_S(e)';
            'EX_3hpp(e)';
            'EX_R_3hpt(e)';
            'EX_S_3hpt(e)';
            'EX_styr(e)';
            'EX_phe_L(e)';
            'EX_phstyr(e)'
           };
end

function subSystems = returnSubSystems()
    subSystems = {' 1,3-Propanediol production';
                  ' 1,4-Butanediol production';
                  ' 2,3-Butanediol production';
                  ' 2,3-Butanediol production';
                  ' 2,3-Butanediol production';
                  ' 3-Hydroxybutyrate production';
                  ' 3-Hydroxybutyrate production';
                  ' 3-Hydroxypropanoate production';
                  ' 3-Hydroxyvalerate production';
                  ' 3-Hydroxyvalerate production';
                  'Styrene production';
                  'Styrene production';
                  'p-Hydroxystyrene production'};
    intermediate = {'Production pathway intermediate synthesis';};
end
