function runOptSwapDed(knockoutNum, swapNum)
% runs OptSwapD
% by Zachary King, 8/13/2012

    disp('running optSwapDed')
    ticID = tic;


    dhRxns = {
        '3OAR100'
        '3OAR40'
        '3OAR60'
        '3OAR80'
        'AGPR'
        'AKGDH'
        'ASAD'
        'DHDPRy'
        'EAR40x'
        'EAR60x'
        'EAR80x'
        'G6PDH2r'
        'GAPD'
        'GLUDy'
        'GND'
        'HPYRRx'
        'HSDy'
        'ICDHyr'
        'IPMD'
        'KARA1'
        'KARA2'
        'MDH'
        'MTHFD'
        'NADH16pp'
        'PDH'
        'PGCD'
        'SHK3Dr'
        'TRSARr'
        'ME1'
        'LDH_D'
        'LCARR'
             };

    % name the run
    run = sprintf('optSwapDed-%ddhs-%dKOs-%dswaps', length(dhRxns), knockoutNum, swapNum);
    logFile = sprintf('%s-%s.csv', run, datestr(now, 'yy-mm-dd_HH_MM_SS'));

    % load the model
    model = loadModelNamed('iJO');

    global biomassRxn
    biomassRxn = model.rxns(model.c~=0);
    global fileId
    fileId = fopen(logFile, 'a');

    myPrint('%s\n', run);


    % set parameters
    targetRxns = {'EX_lac-D(e)'
                  'EX_glyc(e)'
                  'EX_pyr(e)'
                  'EX_fum(e)'
                  'EX_akg(e)'
                 };
    myPrint('dhRxns,', []); printReactions(dhRxns);
    myPrint('\n', []);



    isAerobic = false;
    substrate = 'EX_glc(e)';
    transhydrogenaseKnockout = true;
    % knockoutNum = 0;
    % swapNum = 4;
    minBiomass = 0.1;

    % prepare the model
    if transhydrogenaseKnockout
        % knock out th's
        thRxns = {'NADTRHD', 'THD2pp'};
        model = changeRxnBounds(model, thRxns, 0, 'b');
    end
    model = changeRxnBounds(model, 'EX_glc(e)', 0, 'l');
    model = changeRxnBounds(model, substrate, -20, 'l');
    if isAerobic
        model = changeRxnBounds(model, 'EX_o2(e)', -20, 'l');
    else
        model = changeRxnBounds(model, 'EX_o2(e)', 0, 'l');
        % BAD REACTIONS, DOWN BOY
        model = changeRxnBounds(model, {'CAT';'SPODM';'SPODMpp'}, [0;0;0], 'b');
        % ask Ali about these.
        % can make oxygen...turn off for anaerobic simulations
    end
    model = changeRxnBounds(model, biomassRxn, minBiomass, 'l');

    % save wildtype reaction--with transhydrogenase knockout
    modelWT = model;

    % Swap reactions, keeping wild type dehydrogenases
    [model, newNames] = modelSwap(model, dhRxns, true);

    % Choose reactions for knocking
    ssExcludeList = {'Cell Envelope Biosynthesis',
                     'Exchange',
                     'Inorganic Ion Transport and Metabolism',
                     'Lipopolysaccharide Biosynthesis / Recycling',
                     'Murein Biosynthesis',
                     'Murein Recycling',
                     'Transport, Inner Membrane',
                     'Transport, Outer Membrane',
                     'Transport, Outer Membrane Porin',
                     'tRNA Charging'
                    };
    nCarbonThr = 10;
    [selectedRxns,sets,trimmedSets,geneList,single_KO_GRs,~] = ...
        getOptknockTargets(model,ssExcludeList,nCarbonThr,false);
    % TODO: decide whether to use this REDUCED model
    % TODO: compare essentiality check in getOptknockTargets to RobustKnock/createRobustOneWayFluxesModel

    % set options
    options.useCobraSolver = 0;
    options.knockType = 2; % optSwapD
    options.knockoutNum = knockoutNum;
    options.knockableRxns = selectedRxns;
    options.notKnockableRxns = [dhRxns; newNames];
    options.swapNum = swapNum;
    options.dhRxns = dhRxns;

    myPrint(['target,yield at min biomass,,yield at max biomass,,' ...
             'kos,swaps\n'], []);
    myPrint(',th knockout,optSwap,th knockout,optSwap\n', []);

    for i = 1:length(targetRxns)
        targetRxn = targetRxns{i};

        isSpecial = 1;
        % deal with special cases
        if strcmp(targetRxn, 'EX_h2(e)') % turn on FHL for hydrogen production
            model = changeRxnBounds(model, 'FHL', 1000, 'u');
        elseif strcmp(targetRxn, 'EX_ala-L(e)')  % exporting valine does not matter
            model = changeRxnBounds(model, 'EX_val-L(e)', 0, 'u');
        else
            isSpecial = 0;
        end

        myPrint('%s,', targetRxn);

        disp(sprintf('OptSwapD with target reaction %s', targetRxn));
        options.targetRxn = targetRxn;
        results = optSwapD(model, options);

        knockouts = results.knockoutRxns;
        modelT = modelWT;
        if ~isempty(knockouts)
            modelT = changeRxnBounds(modelT, knockouts, 0, 'b');
        end
        if ~isempty(results.knockoutDhs)
            modelT = modelSwap(modelT, results.knockoutDhs, false);
        end

        printMaxYield(modelWT,targetRxn);
        printMaxYield(modelT, targetRxn);
        printCoupledYield(modelWT, targetRxn);
        printCoupledYield(modelT, targetRxn);

        if ~isempty(knockouts)
            printReactions(knockouts);
        else
            display('no knockouts');
            myPrint(',', []);
        end
        if ~isempty(results.knockoutDhs)
            printReactions(results.knockoutDhs);
        else
            display('no swaps');
            myPrint(',', []);
        end
        myPrint('\n', []);

    end
    t = toc(ticID);
    myPrint('time (min),%.1f', t/60);
    fclose(fileId);
end


function printMaxYield(model, targetRxn)
    model.c = zeros(size(model.c));
    model.c(ismember(model.rxns, targetRxn)) = 1;
    soln = optimizeCbModel(model);
    myPrint('%.4f,', soln.f);
end

function printCoupledYield(model, targetRxn)
    global biomassRxn
    model.c = zeros(size(model.c));
    model.c(ismember(model.rxns, biomassRxn)) = 1;
    soln = optimizeCbModel(model);
    myPrint('%.4f,', soln.x(ismember(model.rxns, targetRxn)));
end

function printReactions(reactions)
    for j=1:length(reactions)
        myPrint('''%s'';', reactions{j});
    end
end

function myPrint(string, val)
    global fileId
    display(sprintf(string, val));
    fprintf(fileId, string, val);
end

% modelT = model;
% modelT.c = zeros(size(model.c));
% modelT.c(ismember(modelT.rxns, targetRxns{i})) = 1;
% modelT.lb(ismember(modelT.rxns, biomassRxn)) = 0.01;
% sol = optimizeCbModel(modelT);
% maxProd{i} = sol.f

% logFile = ['log-' run '.txt'];
% fileId = fopen(logFile, 'a');
% fprintf(fileId, '%s: %g\n', targetRxns{i}, sol.f);
% fclose(fileId);


% % write data file;
% save(sprintf('%s_%d_%s',run, i, datestr(now, 'yy-mm-dd_HH_MM_SS')));
% % write to log file
% t = toc(ticID);
% logFile = ['log-' run '.txt'];
% fileId = fopen(logFile,'a');
% fprintf(fileId,'%s: %s finished in %g min\n',...
%                 datestr(now,'yy-mm-dd_HH_MM_SS'), run, t/60);

% knockoutString = '';
% for j=1:length(results.knockoutRxns)
%     knockoutString = [knockoutString ' ' results.knockoutRxns{j}];
% end
% display(['knockouts: ', knockoutString]);
% fprintf(fileId,'knocked rxns: %s\n', knockoutString);


% newModel = changeRxnBounds(model, results.knockoutRxns, 0, 'b');
% compareOptions.substrateRxn = 'EX_glc(e)';
% compareOptions.targetRxn = options.targetRxn;
% compareOptions.isAerobic = isAerobic;
% compareOptions.pairwiseCompare = true;
% compareOptions.showEnvelope = true;
% output = compareModels({newModel, model}, compareOptions)
% % TODO: print to file
% % fprintf(fileId,'%s\n',output);

% fclose(fileId);



% try
%     tw.updateStatus(sprintf('@zakandrewking %s finished', run));
% end
