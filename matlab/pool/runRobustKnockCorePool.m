function runRobustKnockCorePool(opt)
% runs OptSwap
% by Zachary King, 8/13/2012

    disp('running OptKnock')
    ticID = tic; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set default parameters
    if ~exist('opt','var')
        opt = struct();
    end
    if ~isfield(opt,'knockoutNum'), opt.knockoutNum = 4; end
    if ~isfield(opt, 'targetRxn') opt.targetRxn = 'EX_succ(e)'; end
    if ~isfield(opt, 'startWithKnocks'), opt.startWithKnocks = []; end
    if ~isfield(opt, 'experiment'), opt.experiment = 'succinate_core'; end
    if ~isfield(opt, 'notes'), opt.notes = ''; end
    if ~isfield(opt, 'logFile'), opt.logFile = 'core_optknock_pool_log.csv'; end
    if ~isfield(opt, 'aerobicString'), opt.aerobicString = 'anaerobic'; end
    if ~isfield(opt, 'solverParams'), opt.solverParams = struct(); end %seconds
    if ~isfield(opt, 'substrate'), opt.substrate = 'EX_glc(e)'; end
    if ~isfield(opt, 'useCobraSolver'), opt.useCobraSolver = true; end
    
    % make variables local
    knockoutNum = opt.knockoutNum;
    targetRxn = opt.targetRxn;
    startWithKnocks = opt.startWithKnocks;
    experiment = opt.experiment;
    notes = opt.notes;
    logFile = opt.logFile;
    aerobicString = opt.aerobicString;
    solverParams = opt.solverParams;
    substrate = opt.substrate;
    useCobraSolver = opt.useCobraSolver;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global biomassRxn minBiomass
    minBiomass = 0.05;
    [model, biomassRxn] = setupModel('ecoli_core_model', substrate, aerobicString);
    model = changeRxnBounds(model, substrate, -20, 'l');

    % edit model with starting swaps/knocks
    if ~isempty(startWithKnocks)
        model = changeRxnBounds(model, startWithKnocks, 0, 'b');
    end

    global fileId
    fileId = fopen(logFile, 'a');
    
    % set options
    options.numDel = knockoutNum;
    options.targetRxn = targetRxn;
    constrOpt.rxnList = [biomassRxn, {'ATPM'}];
    constrOpt.values = [minBiomass, 8.39];
    constrOpt.sense = 'GE';
    % ignore exchange, transport, and ATPM and biomass
    selectedRxnList = model.rxns(~findExcRxns(model) & ...
                                 ~ismember(model.rxns, constrOpt.rxnList) & ...
                                 cellfun(@isempty, strfind(model.rxnNames, 'transport')));

    prevSolutions = {};
    maxFail = 20; fails = 0;
    while true
        myPrint('\n',[]);
        lTic = tic;

        modelTR = model;
        modelTR = setupModelForTarget(modelTR, targetRxn);

        fprintf('OptKnock with target reaction %s\n', targetRxn);
        error('Unfinished');
        optKnockSol = OptSwap(modelTR, selectedRxnList, options, constrOpt, prevSolutions);

        if ~strcmp(optKnockSol.origStat, 'OPTIMAL') 
            if strcmp(optKnockSol.origStat, 'INTERRUPTED') 
                break
            end
            myPrint('%s,', optKnockSol.origStat);
            if fails >= maxFail
                break
            else
                fails = fails + 1;
                continue
            end
        end
        
        knockouts = optKnockSol.rxnList;

        modelWT = setupModelForTarget(model, targetRxn);
        modelT = modelWT;
        if ~isempty(knockouts)
            modelT = changeRxnBounds(modelT, knockouts, 0, 'b');
        end

        t = toc(lTic);

        % Project   Date    Time    Name    Script  Knockout num    Swap num
        % Intervention num ... 
        myPrint('OptKnock,', []);
        myPrint('%s,', experiment);
        myPrint('%s,', datestr(now, 'mm/dd/yyyy'));
        myPrint('%s,', datestr(now, 'HH:MM:SS'));
        myPrint('optSwap.m,', []);
        myPrint('%d,', knockoutNum);
        myPrint('%s,', targetRxn);
        printMaxYield(modelWT,targetRxn);
        printMaxYield(modelT, targetRxn);
        % solnWT = printCoupledYield(modelWT, targetRxn);
        solnT = printCoupledYield(modelT, targetRxn, 1);        
        solnT = printCoupledYield(modelT, targetRxn, -1);
        if ~isempty(knockouts)
            printReactions(knockouts);
        else
            display('no knockouts');
            myPrint(',', []);
        end
        myPrint('%.2g,', t);
        myPrint('%.1f,', minBiomass);
        myPrint('%s,', aerobicString); % Aerobicity	
        myPrint('%s,', substrate);     % Substrate
        myPrint('%s,', optKnockSol.origStat);

        if isempty(knockouts) 
            if fails >= maxFail
                break
            else
                fails = fails + 1;
            end
        end

        prevSolutions{end+1} = knockouts;

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

function soln = printCoupledYield(model, targetRxn, sense)
    global biomassRxn
    model = changeObjective(model, targetRxn, sense);
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
