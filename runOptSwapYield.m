function runOptSwapYield(options)

    substrate = options.substrate;
    aerobicString = options.aerobicString;
    swapNum = options.swapNum;
    logFile = options.logFile;
    if isfield(options, 'dhRxns')
        dhRxns = options.dhRxns
    else
        dhRxns = dhRxnList(30);
    end

    global fileId
    fileId = fopen(logFile, 'a');
    fprintf(fileId, 'target\taerobic\tsubstrate\tthko\tswaps\n')

    targetRxns = returnTargetRxns();    
    [model,biomassRxn] = setupModel('iJO',substrate,aerobicString,'thko');
    
    for i=1:length(targetRxns)
        modelT = model;
        opt.targetRxn = targetRxns{i};
        modelT = setupModelForTarget(modelT, opt.targetRxn);
        opt.useCobraSolver = true;
        opt.biomassRxn = biomassRxn;
        opt.swapNum = swapNum;
        opt.minBiomass = 0.1;
        opt.dhRxns = dhRxns;
        results = optSwapYield(modelT, opt);
        myPrint('%s\t',targetRxns{i});
        myPrint('%s\t', aerobicString);
        myPrint('%s\t', substrate);
        myPrint('%d\t', 1);
        myPrint('%.4f\t', results.f_k);
        if ~isempty(results.knockoutDhs)
            printReactions(results.knockoutDhs);
        else
            display('no swaps');
            myPrint('\t', []);
        end
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
             'DM_5DRIB';
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
             'DM_OXAM';
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
             'DM_MTHTHF';
             'EX_fum(e)';
             'EX_quin(e)';
             'EX_dha(e)';
             'EX_h(e)';
             'EX_co2(e)';
             'EX_12ppd-R(e)';
             'EX_succ(e)';
             'EX_for(e)';};
end