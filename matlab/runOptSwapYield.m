function runOptSwapYield(options)

    substrate = options.substrate;
    aerobicString = options.aerobicString;
    swapNum = options.swapNum;
    logFile = options.logFile;
    thko = options.thko;
    if isfield(options, 'dhRxns')
        dhRxns = options.dhRxns;
    else
        dhRxns = dhRxnList('yield');
    end
    if isfield(options, 'targetRxns')
        targetRxns = options.targetRxns;
    else
        targetRxns = returnTargetRxns();    
    end        
    if isfield(options, 'modelname')
        modelname = options.modelname;
    else
        modelname = 'iJO';
    end        
    if isfield(options, 'minBiomass')
        minBiomass = options.minBiomass;
    else
        minBiomass = 0.1
    end
    if isfield(options, 'nondefaultReactionBounds')
        nondefaultReactionBounds = options.nondefaultReactionBounds;
    else
        nondefaultReactionBounds = [];
    end

    global fileId
    fileId = fopen(logFile, 'a');

    [model,biomassRxn] = setupModel(modelname,substrate,aerobicString,thko);
    for i=1:length(targetRxns)
        lTic = tic;
        opt.useCobraSolver = true;
        opt.biomassRxn = biomassRxn;
        opt.swapNum = swapNum;
        opt.minBiomass = minBiomass;
        opt.dhRxns = dhRxns;
        opt.targetRxn = targetRxns{i};
        
        modelT = model;
        modelT = setupModelForTarget(modelT, opt.targetRxn);
        for z=1:length(nondefaultReactionBounds)
            modelT = changeRxnBounds(modelT, nondefaultReactionBounds{z}{1}, ...
                                     nondefaultReactionBounds{z}{2}, 'l');
            modelT = changeRxnBounds(modelT, nondefaultReactionBounds{z}{1}, ...
                                     nondefaultReactionBounds{z}{3}, 'u');
        end

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
        modelT = changeRxnBounds(modelT, biomassRxn, minBiomass, 'l');
        soln = optimizeCbModel(modelT);
        
        myPrint('%s\t',targetRxns{i});
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
    rxns =  {'EX_h2s(e)';
             'EX_cys-L(e)';
             'EX_hom-L(e)';
             'EX_thr-L(e)';
             'EX_15dap(e)';
             'EX_lys-L(e)';
             'EX_ile-L(e)';
             'EX_leu-L(e)';
             'EX_gly(e)';
             'EX_5mtr(e)';
             'EX_spmd(e)';
             'EX_pro-L(e)';
             'EX_cgly(e)';
             'EX_asp-L(e)';
             'EX_gthrd(e)';
             'EX_ser-L(e)';
             'EX_thym(e)';
             'EX_ptrc(e)';
             'EX_glyc(e)';
             'EX_pyr(e)';
             'EX_orn(e)';
             'EX_agm(e)';
             'EX_urea(e)';
             'EX_arg-L(e)';
             'EX_glyald(e)';
             'EX_lac-L(e)';
             'EX_h2o(e)';
             'EX_asn-L(e)';
             'EX_thymd(e)';
             'EX_LalaDgluMdapDala(e)';
             'EX_acser(e)';
             'EX_ac(e)';
             'EX_LalaDgluMdap(e)';
             'EX_phe-L(e)';
             'EX_tyr-L(e)';
             'EX_ala-D(e)';
             'EX_alaala(e)';
             'EX_12ppd-S(e)';
             'EX_his-L(e)';
             'EX_g3pe(e)';
             'EX_g3pg(e)';
             'EX_enter(e)';
             'EX_feenter(e)';
             'EX_4abut(e)';
             'EX_ura(e)';
             'EX_ade(e)';
             'EX_pheme(e)';
             'EX_adn(e)';
             'EX_glyc3p(e)';
             'EX_etha(e)';
             'EX_glu-L(e)';
             'EX_uri(e)';
             'EX_indole(e)';
             'EX_trp-L(e)';
             'EX_cytd(e)';
             'EX_gua(e)';
             'EX_hxan(e)';
             'EX_xan(e)';
             'EX_alltn(e)';
             'EX_acolipa(e)';
             'EX_xtsn(e)';
             'EX_colipa(e)';
             'EX_ins(e)';
             'EX_kdo2lipid4(e)';
             'EX_acald(e)';
             'EX_enlipa(e)';
             'EX_anhgm(e)';
             'EX_eca4colipa(e)';
             'EX_lipa_cold(e)';
             'EX_lipa(e)';
             'EX_colipap(e)';
             'EX_glyc-R(e)';
             'EX_cit(e)';
             'EX_akg(e)';
             'EX_h2(e)';
             'EX_glyclt(e)';
             'EX_mal-L(e)';
             'EX_fe2(e)';
             'EX_fum(e)';
             'EX_quin(e)';
             'EX_dha(e)';
             'EX_h(e)';
             'EX_co2(e)';
             'EX_12ppd-R(e)';
             'EX_succ(e)';
             'EX_for(e)';};
end