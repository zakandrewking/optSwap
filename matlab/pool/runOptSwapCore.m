function runOptSwapCore(opt)
% runs OptSwap
% by Zachary King, 8/13/2012

    disp('running optSwap')
    ticID = tic; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set default parameters
    if ~exist('opt','var')
        opt = struct();
    end
    if ~isfield(opt,'knockoutNum'), opt.knockoutNum = 5; end
    if ~isfield(opt, 'swapNum'), opt.swapNum = 0; end
    if ~isfield(opt, 'interventionNum'), opt.interventionNum = -1; end
    if ~isfield(opt, 'targetRxns')
        opt.targetRxns = {'EX_pyr(e)'};
    end
    if ~isfield(opt, 'startWithKnocks'), opt.startWithKnocks = []; end
    if ~isfield(opt, 'startWithSwaps'), opt.startWithSwaps = []; end
    if ~isfield(opt, 'experiment'), opt.experiment = 'succinate_core'; end
    if ~isfield(opt, 'notes'), opt.notes = ''; end
    if ~isfield(opt, 'logFile'), opt.logFile = 'core_robustknock_log.csv'; end
    if ~isfield(opt, 'swapAllDhs'), opt.swapAllDhs = false; end
    if ~isfield(opt, 'rxnSet'), opt.rxnSet = {}; end
    if ~isfield(opt, 'aerobicString'), opt.aerobicString = 'anaerobic'; end
    if ~isfield(opt, 'solverParams'), opt.solverParams = struct(); end %seconds
    if ~isfield(opt, 'substrate'), opt.substrate = 'EX_glc(e)'; end
    if ~isfield(opt, 'printIntermediateSolutions') % unfinished
        opt.printIntermediateSolutions = false; 
    end
    if ~isfield(opt, 'useCobraSolver'), opt.useCobraSolver = true; end
    if ~isfield(opt, 'canKnockDHs'), opt.canKnockDHs = false; end
    if ~isfield(opt, 'knockType'), opt.knockType = 0; end
    if ~isfield(opt, 'allowDehydrogenaseKnockout')
        opt.allowDehydrogenaseKnockout = true; 
    end
    
    % make variables local
    knockoutNum = opt.knockoutNum;
    swapNum = opt.swapNum;
    interventionNum = opt.interventionNum;
    targetRxns = opt.targetRxns;
    startWithKnocks = opt.startWithKnocks;
    startWithSwaps = opt.startWithSwaps;
    experiment = opt.experiment;
    notes = opt.notes;
    logFile = opt.logFile;
    swapAllDhs = opt.swapAllDhs;
    rxnSet = opt.rxnSet;
    aerobicString = opt.aerobicString;
    solverParams = opt.solverParams;
    substrate = opt.substrate;
    printIntermediateSolutions = opt.printIntermediateSolutions;
    useCobraSolver = opt.useCobraSolver;
    canKnockDHs = opt.canKnockDHs;
    knockType = opt.knockType;
    allowDehydrogenaseKnockout = opt.allowDehydrogenaseKnockout;

    % check values
    if ~iscell(targetRxns)
        targetRxns = {targetRxns};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % allow knockout of dhRxns
    if canKnockDHs
        if swapNum>0, error('cannot knock dehydrogenases if swapsNum > 0'); end
        dhRxns = [];
        knockableRxns = dhRxnList(21);
        fprintf('RoubstKnock--Can knockout dehydrogenases\n');
    else
        % dhRxns = dhRxnList(21); 
        dhRxns = [];
        knockableRxns = {}; 
    end
    
    fprintf('%d dehydrogenase reactions\n', length(dhRxns));

    % name the run
    run = sprintf('%ddhs-%dKOs-%dswaps', length(dhRxns),...
                  knockoutNum, swapNum);

    global biomassRxn minBiomass
    minBiomass = 0.05;
    [model, biomassRxn] = setupModel('ecoli_core_model', substrate, aerobicString);
  
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
    
    % set options
    options.knockType = knockType;
    options.knockoutNum = knockoutNum;
    options.interventionNum = interventionNum;
    options.knockableRxns = model.rxns(~findExcRxns(model));
    options.notKnockableRxns = model.rxns(findExcRxns(model));
    options.swapNum = swapNum;
    options.dhRxns = dhRxns;
    options.solverParams = solverParams;
    options.useCobraSolver = useCobraSolver;
    options.allowDehydrogenaseKnockout = allowDehydrogenaseKnockout;

    for i = 1:length(targetRxns)
        myPrint('\n',[]);
        lTic = tic;
        targetRxn = targetRxns{i};

        modelTR = model;
        modelTR = setupModelForTarget(modelTR, targetRxn);

        fprintf('OptSwap with target reaction %s\n', targetRxn);
        options.targetRxn = targetRxn;
        options
        results = optSwap(modelTR, options);

        knockouts = results.knockoutRxns;

        modelWT = setupModelForTarget(model, targetRxn);
        modelT = modelWT;
        if ~isempty(knockouts)
            modelT = changeRxnBounds(modelT, knockouts, 0, 'b');
        end
        if ~isempty(results.swapRxns)
            modelT = modelSwap(modelT, results.swapRxns, false);
        end

        t = toc(lTic);

        % Project   Date    Time    Name    Script  Knockout num    Swap num
        % Intervention num ... 
        myPrint('OptSwap,', []);
        myPrint('%s,', experiment);
        myPrint('%s,', datestr(now, 'mm/dd/yyyy'));
        myPrint('%s,', datestr(now, 'HH:MM:SS'));
        myPrint('%s,', run);
        myPrint('optSwap.m,', []);
        myPrint('%d,', options.knockoutNum);
        myPrint('%d,', options.swapNum);
        myPrint('%d,', options.interventionNum); 
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
        if ~isempty(results.swapRxns)
            printReactions(results.swapRxns);
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
        if isnumeric(results.inform)
            myPrint('%d,', results.inform);
        else
            myPrint('%s,', results.inform);
        end
        myPrint('%s,', results.solver);
        fi = fieldnames(solverParams);
        for i=1:length(fi)
            myPrint('%s: ', fi{i});
            myPrint('%.2e,', getfield(solverParams, fi{i}));
        end
    end
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
