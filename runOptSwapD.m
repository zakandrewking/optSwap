function runOptSwapD(opt)
% runs OptSwapD
% by Zachary King, 8/13/2012



    disp('running optSwapD')
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
    
    % make variables local
    knockoutNum = opt.knockoutNum;
    swapNum = opt.swapNum;
    targetRxns = opt.targetRxns;
    startWithKnocks = opt.startWithKnocks;
    startWithSwaps = opt.startWithSwaps;
    experiment = opt.experiment;
    notes = opt.notes;
    
    % check values
    if ~iscell(targetRxns)
        targetRxns = {targetRxns};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
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
    run = sprintf('%ddhs-%dKOs-%dswaps', length(dhRxns),...
                  knockoutNum, swapNum);
    logFile = 'database.csv';

    global biomassRxn minBiomass
    [model, biomassRxn] = setupModel('iJO','EX_glc(e)',false,true);
    minBiomass = 0.1;
    
    % edit model with starting swaps/knocks
    if ~isempty(startWithKnocks)
        model = changeRxnBounds(model, startWithKnocks, 0, 'b');
    end
    if ~isempty(startWithSwaps)
        model = modelSwap(model, startWithSwaps, false);
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
    reducedModelFilename = sprintf('reducedModel-%s-%s.mat', 'iJO', 'glc');
    % need a new reduced model if we start with knocks or swaps
    noStartProc = isempty(startWithKnocks) && isempty(startWithSwaps);
    if exist(reducedModelFilename,'file') == 2 && noStartProc 
        load(reducedModelFilename);
    else
        [selectedRxns,~,~,~,~,reducedModel] = ...
            getOptknockTargets(model,ssExcludeList,nCarbonThr,true);
        if noStartProc
            save(reducedModelFilename,'selectedRxns','reducedModel');
        end
    end
    
    % remove dhRxns from set if they are not in reduced model
    dhRxns(~ismember(dhRxns, reducedModel.rxns)) = [];
    
    % set options
    options.useCobraSolver = 0;
    options.knockType = 2; % optSwapD
    options.knockoutNum = knockoutNum;
    % options.knockableRxns = {};
    notSelectedRxns = model.rxns(~ismember(model.rxns,selectedRxns));
    options.notKnockableRxns = notSelectedRxns;
    options.swapNum = swapNum;
    options.dhRxns = dhRxns;

    for i = 1:length(targetRxns)
        lTic = tic;
        targetRxn = targetRxns{i};

        modelTR = reducedModel; 
        modelTR = setupModelForTarget(modelTR, targetRxn);
        
        fprintf('OptSwap with target reaction %s\n', targetRxn);
        options.targetRxn = targetRxn;
        results = optSwapD(modelTR, options);

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

        % Project	Date	Time	Name	Script	Knockout num	Swap num	
        % Target reaction	Yield at min biomass / TH KO	Yield at min biomass
        %  / OptSwap	Yield at max biomass  / TH KO	Yield at max biomass / 
        % OptSwap,Knockouts,Swaps,Time (min),Min biomass,Dehydrogenase reactions
        myPrint('OptSwap,', []);
        myPrint('%s,', experiment);
        myPrint('%s,', datestr(now, 'mm/dd/yyyy'));
        myPrint('%s,', datestr(now, 'HH:MM:SS'));
        myPrint('%s,', run);
        myPrint('optSwapD.m,', []);
        myPrint('%d,', options.knockoutNum);
        myPrint('%d,', options.swapNum);
        myPrint('%s,', targetRxn);
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
        printReactions(dhRxns); 
        myPrint('%s,', notes); %notes
        myPrint('%d,', results.exitFlag);
        myPrint('%d,', results.inform);
        myPrint('\n', []);

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
