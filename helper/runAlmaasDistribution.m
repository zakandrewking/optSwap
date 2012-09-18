function runAlmaasDistribution

    tw = setupTwitty;
    model = loadModelNamed('iAF1260b');
    options.subs = 'EX_glc(e)';
    % options.possibleLoopRxns = {'TRSARr','HPYRRx'};
    options.autRemLoops = true;
    % options.usePFBA = true;
    [returnRxns,fluxes] = almaasDistribution(model,options);
    print = [returnRxns, num2cell(fluxes)]
    % save_to_base(1);
    S = tw.updateStatus('this gremlin is awake');
end