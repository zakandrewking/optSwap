function runOptSwap(opt)
% runs OptSwap
% by Zachary King, 8/13/2012

    disp('running optSwap')
    ticID = tic; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set default parameters
    if ~exist('opt','var')
        opt = struct();
    end
    if ~isfield(opt,'knockoutNum'), opt.knockoutNum = 0; end
    if ~isfield(opt, 'swapNum'), opt.swapNum = 0; end
    if ~isfield(opt, 'targetRxns')
        opt.targetRxns = {
            'EX_etoh(e)',
            'EX_for(e)',
            'EX_succ(e)',
            'EX_ac(e)',
                         };
    end
    if ~isfield(opt, 'startWithKnocks'), opt.startWithKnocks = []; end
    if ~isfield(opt, 'startWithSwaps'), opt.startWithSwaps = []; end
    if ~isfield(opt, 'experiment'), opt.experiment = 'no name'; end
    if ~isfield(opt, 'notes'), opt.notes = ''; end
    if ~isfield(opt, 'logFile'), opt.logFile = 'database.csv'; end
    if ~isfield(opt, 'swapAllDhs'), opt.swapAllDhs = false; end
    if ~isfield(opt, 'rxnSet'), opt.rxnSet = {}; end
    if ~isfield(opt, 'aerobicString'), opt.aerobicString = 'anaerobic'; end
    if ~isfield(opt, 'maxTime'), opt.maxTime = 60; end
    if ~isfield(opt, 'substrate'), opt.substrate = 'EX_glc(e)'; end
    if ~isfield(opt, 'printIntermediateSolutions')
        opt.printIntermediateSolutions = false; 
    end
    
    % make variables local
    knockoutNum = opt.knockoutNum;
    swapNum = opt.swapNum;
    targetRxns = opt.targetRxns;
    startWithKnocks = opt.startWithKnocks;
    startWithSwaps = opt.startWithSwaps;
    experiment = opt.experiment;
    notes = opt.notes;
    logFile = opt.logFile;
    swapAllDhs = opt.swapAllDhs;
    rxnSet = opt.rxnSet;
    aerobicString = opt.aerobicString;
    maxTime = opt.maxTime;
    substrate = opt.substrate;
    printIntermediateSolutions = opt.printIntermediateSolutions;

    % check values
    if ~iscell(targetRxns)
        targetRxns = {targetRxns};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    dhRxns = dhRxnList(31);

    % name the run
    run = sprintf('%ddhs-%dKOs-%dswaps', length(dhRxns),...
                  knockoutNum, swapNum);

    global biomassRxn minBiomass
    minBiomass = 0.1;
    [model, biomassRxn] = setupModel('iJO',substrate,aerobicString,'thko');
  
    % edit model with starting swaps/knocks
    if ~isempty(startWithKnocks)
        model = changeRxnBounds(model, startWithKnocks, 0, 'b');
    end
    if ~isempty(startWithSwaps)
        model = modelSwap(model, startWithSwaps, false);
    end
    if swapAllDhs
        model = modelSwap(model, dhRxns, true);
    end

    global fileId
    fileId = fopen(logFile, 'a');

    % Choose reactions for knocking
    ssExcludeList = {'Cell Envelope Biosynthesis',...
                     'Exchange',...
                     'Inorganic Ion Transport and Metabolism',...
                     'Lipopolysaccharide Biosynthesis / Recycling',...
                     'Murein Biosynthesis',...
                     'Murein Recycling',...
                     'Transport, Inner Membrane',...
                     'Transport, Outer Membrane',...
                     'Transport, Outer Membrane Porin',...
                     'tRNA Charging'
                    };
    nCarbonThr = 10;

    % load or make reduced model
    if isempty(rxnSet)
        reducedModelFilename = sprintf('reducedModel-%s-%s.mat', 'iJO', ...
                                       substrate);
    else
        reducedModelFilename = sprintf('reducedModel-%s-%s-%d.mat', 'iJO', ...
                                       substrate, length(rxnSet));
    end
    % need a new reduced model if we start with knocks or swaps
    noStartProc = isempty(startWithKnocks) && isempty(startWithSwaps);
    if exist(reducedModelFilename,'file') == 2 && noStartProc
        load(reducedModelFilename);
    else
        if isempty(rxnSet)
            [selectedRxns,~,~,~,~,reducedModel] = ...
                getOptknockTargets(model,ssExcludeList,nCarbonThr,true);
        else
            selectedRxns = rxnSet;
            reducedModel = reduceModel(model);
        end
        if noStartProc
            save(reducedModelFilename,'selectedRxns','reducedModel');
        end
    end
    fprintf('\nNumber of selected reactions: %d\n', length(selectedRxns))

    % remove dhRxns from set if they are not in reduced model
    dhRxns(~ismember(dhRxns, reducedModel.rxns)) = [];

    % set options
    options.useCobraSolver = 0;
    options.knockType = 2; % optSwap
    options.knockoutNum = knockoutNum;
    % options.knockableRxns = {};
    notSelectedRxns = reducedModel.rxns(~ismember(reducedModel.rxns,selectedRxns));
    options.notKnockableRxns = notSelectedRxns;
    options.swapNum = swapNum;
    options.dhRxns = dhRxns;
    options.maxTime = maxTime;
    options.printIntermediateSolutions = printIntermediateSolutions;
    if printIntermediateSolutions
        options.intermediateSolutionsFile = [experiment '-MILPsols.csv'];
    end

    for i = 1:length(targetRxns)
        myPrint('\n',[]);
        lTic = tic;
        targetRxn = targetRxns{i};

        modelTR = reducedModel;
        modelTR = setupModelForTarget(modelTR, targetRxn);

        fprintf('OptSwap with target reaction %s\n', targetRxn);
        options.targetRxn = targetRxn;
        results = optSwap(modelTR, options);

        knockouts = results.knockoutRxns;

        modelWT = setupModelForTarget(model, targetRxn);
        modelT = modelWT;
        if ~isempty(knockouts)
            modelT = changeRxnBounds(modelT, knockouts, 0, 'b');
        end
        if ~isempty(results.knockoutDhs)
            modelT = modelSwap(modelT, results.knockoutDhs, false);
        end

        t = toc(lTic);

        % Project   Date    Time    Name    Script  Knockout num    Swap num
        % Target reaction   Yield at min biomass / TH KO    Yield at min biomass
        %  / OptSwap    Yield at max biomass  / TH KO   Yield at max biomass /
        % OptSwap,Knockouts,Swaps,Time (min),Min biomass,Dehydrogenase reactions
        myPrint('OptSwap,', []);
        myPrint('%s,', experiment);
        myPrint('%s,', datestr(now, 'mm/dd/yyyy'));
        myPrint('%s,', datestr(now, 'HH:MM:SS'));
        myPrint('%s,', run);
        myPrint('optSwap.m,', []);
        myPrint('%d,', options.knockoutNum);
        myPrint('%d,', options.swapNum);
        myPrint('%s,', targetRxn);
        myPrint('%.4f,', results.f_k);
        printMaxYield(modelWT,targetRxn);
        printMaxYield(modelT, targetRxn);
        solnWT = printCoupledYield(modelWT, targetRxn);
        solnT = printCoupledYield(modelT, targetRxn);
        printSsp(modelWT, targetRxn, solnWT);
        printSsp(modelT, targetRxn, solnT);
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
        myPrint('%.1f,', t/60);
        myPrint('%.1f,',minBiomass);
        myPrint('%s,', aerobicString); % Aerobicity	
        myPrint('%s,', substrate);     % Substrate
        printReactions(dhRxns);
        if ~isempty(startWithSwaps) || ~isempty(startWithKnocks)
            printReactions([startWithSwaps;startWithKnocks]);
        else
            myPrint(',', []);
        end
        myPrint('%d,', results.exitFlag);
        myPrint('%d,', results.inform);
    end
    % t = toc(ticID);
    % myPrint('time (min),%.1f', t/60);
    fclose(fileId);
end


function printMaxYield(model, targetRxn)
    global biomassRxn minBiomass
    model = changeRxnBounds(model, biomassRxn, minBiomass, 'l');
    model.c = zeros(size(model.c));
    model.c(ismember(model.rxns, targetRxn)) = 1;
    soln = optimizeCbModel(model);
    myPrint('%.4f,', soln.f);
end

function soln = printCoupledYield(model, targetRxn)
    global biomassRxn
    model.c = zeros(size(model.c));
    model.c(ismember(model.rxns, biomassRxn)) = 1;
    soln = optimizeCbModel(model);
    if ~isempty(soln.x)
        myPrint('%.4f,', soln.x(ismember(model.rxns, targetRxn)));
    else
        myPrint('no sol,', []);
    end
end

function printSsp(model, targetRxn, prevSoln)
    if nargin < 3
        global biomassRxn
        model.c = zeros(size(model.c));
        model.c(ismember(model.rxns, biomassRxn)) = 1;
        soln = optimizeCbModel(model);
    else
        soln = prevSoln;
    end
    if ~isempty(soln.x)
        myPrint('%.4f,', soln.x(ismember(model.rxns, targetRxn))*soln.f);
    else
        myPrint('no sol,', []);
    end
end

function printReactions(reactions)
    for j=1:length(reactions)
        myPrint('''%s'';', reactions{j});
    end
    myPrint(',', []);
end

function myPrint(string, val)
    global fileId
    display(sprintf(string, val));
    fprintf(fileId, string, val);
end
