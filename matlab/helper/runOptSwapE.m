function runOptSwapE


    setupCobraSolver;
    tw = setupTwitty;


    % run optSwapE
    model = loadModelNamed('iAF1260b');
    thRxns = {'NADTRHD', 'THD2pp'};

    % options.dhRxns = {'TRDR'};
    options.substrateExRxns = {'EX_glc(e)'};
    options.isAer = 0;
    options.targetRxns = {'EX_ac(e)', 'EX_for(e)', 'EX_lac_D(e)',...
                        'EX_etoh(e)', 'EX_succ(e)' };
    options.thKO = true; 
    options.filename = 'mixed-ferm-ThKO-doubles';
    optSwapE(model, thRxns, options)

    % save_to_base(1);
    tw.updateStatus(sprintf('@zakandrewking %s finished', options.filename)); 
end