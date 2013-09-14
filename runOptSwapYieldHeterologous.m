function runOptSwapYieldHeterologous(options)

    substrate = options.substrate;
    aerobicString = options.aerobicString;
    swapNum = options.swapNum;
    logFile = options.logFile;
    thko = options.thko;
    if isfield(options, 'targetRxns')
        targetRxns = options.targetRxns;
    else
        targetRxns = returnTargetRxns();    
    end        

    global fileId
    fileId = fopen(logFile, 'a');
    
    dhRxns = dhRxnList('yield');

    subSystems = returnSubSystems();
    [model,biomassRxn] = setupModel('iJO-h',substrate,aerobicString,thko);
    for i=1:length(targetRxns)
        lTic = tic;
        opt.useCobraSolver = true;
        opt.biomassRxn = biomassRxn;
        opt.swapNum = swapNum;
        opt.minBiomass = 0.1;
        opt.dhRxns = dhRxns;
        opt.targetRxn = targetRxns{i};

        modelT = model;
        modelT = setupModelForTarget(modelT, opt.targetRxn);
        modelT = turnOnSubSystem(modelT, i, subSystems);
        
        if opt.swapNum==0
            f_k = '';
            knockoutDhs = {};
        else
            results = optSwapYield(modelT, opt);
            f_k = sprintf('%.4f', results.f_k);
            knockoutDhs = results.knockoutDhs;
        end
        
        modelT = modelSwap(modelT, knockoutDhs, false);
        modelT = changeObjective(modelT, opt.targetRxn);
        modelT = changeRxnBounds(modelT, biomassRxn, 0.1, 'l');
        soln = optimizeCbModel(modelT);
       
        myPrint('%s\t',targetRxns{i});
        myPrint('%s\t', subSystems{i});
        myPrint('%s\t', aerobicString);
        myPrint('%s\t', substrate);
        myPrint('%d\t', swapNum);
        myPrint('%s\t', thko);
        myPrint('%s\t', f_k);
        myPrint('%.4f\t', soln.f);
        if ~isempty(knockoutDhs)
            printReactions(knockoutDhs);
        else
            display('no swaps');
            myPrint('\t', []);
        end
        t = toc(lTic);
        myPrint('%.1f', t);
        myPrint('\n',[]);
        % add swaps list
    end
    fclose(fileId);
end

function myPrint(string, val)
    global fileId
    fprintf(string, val);
    fprintf(fileId, string, val);
end

function printReactions(reactions)
    for j=1:length(reactions)
        myPrint('''%s'';', reactions{j});
    end
    myPrint('\t', []);
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
