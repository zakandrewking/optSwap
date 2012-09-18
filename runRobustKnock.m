function runRobustKnock

% chemicalInd - index of chemical to produce
% objectiveInd - organizms objective
% knockoutNum - number of knockouts
% maxW - maximal value of dual variables (higher number will be more
% accurate but takes more calculation time)
% optKnockFlag - indicates if optKnock will be calculated 
% robustKnockFlag - indicates if optKnock will be calculated 
% fixedGlu - size of fixed glucose reaction

    model = loadModelNamed('iAF1260b');        

    almaasOptions = struct('subs','EX_glc(e)',...
                           'possibleLoopRxns',{{'TRSARr','HPYRRx'}});
    dhRxns = almaasDistribution(model,almaasOptions);

    model = modelSwap(model,dhRxns,1);
    
    targetRxn = 'EX_for(e)';
    
    chemicalInd = find(ismember(model.rxns, targetRxn));
    objectiveInd = find(model.c==1);
    knockoutNum = 3;
    maxW = 1e7;
    optKnockFlag = 0;
    robustKnockFlag = 1;
    fixedGlu = 10;
    robustKnock(chemicalInd, objectiveInd, knockoutNum, maxW,...
                optKnockFlag, robustKnockFlag, fixedGlu);
    
end