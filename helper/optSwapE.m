function allSolns = optSwapE(model,thRxns,options)

% OPTSWAPE
%
% Swaps electron carrier specificity for dehydrogenase reactions. Runs FBA
% on model and monitors target rxn outputs (e.g. mixed fermentation products)
%
% Starts with set of rxn knockouts.
% Knocks out transhydrogenases.
%    Knocks out a dehydrogenase rxn and replaces with swapped rxn.
%    Runs productions envelope for ethanol production.
%    Repeats for each dehydrogenase rxn.
%
%
% INPUTS
% model - Model of reconstruction.
% thRxns - transhydrogenase reacions
% options
%   dhRxns - List of dehydrogenase reactions to swap.
%   exUptakeMax - default 20
%   O2Uptake - default 20
%   filename
%
%   isAer - 0 anaerobic, 1 aerobic (default), 2 try both
%   thKO - 0 don't knock out TH rxns (default), 1 knock out TH rxns
%   substrateExRxns - Available carbon sources. default is glucose.
%   targetRxns - default is etoh
%   knockoutRxns - default is none
%
%
% OUTPUT
% solns
%
% Zachary King 7/2/2012
%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check options & set defaults


    if ~isfield(options,'dhRxns')
        almaasOptions = struct('subs','EX_glc(e)',...
                               'possibleLoopRxns',{{'TRSARr','HPYRRx'}});
        options.dhRxns = ...
            almaasDistribution(model,almaasOptions);
    end
    if ~isfield(options,'exUptakeMax')
        options.exUptakeMax = 20;
    end
    if ~isfield(options,'O2Uptake')
        options.O2Uptake = 20;
    end
    if ~isfield(options,'filename')
        options.filename = '';
    end

    if ~isfield(options,'targetRxns')
        options.targetRxns = {'EX_etoh(e)'};
    elseif ~iscell(options.targetRxns)
        options.targetRxns = {options.targetRxns};
    end

    if ~isfield(options,'substrateExRxns')
        options.substrateExRxns = {'EX_glc(e)'};
    elseif ~iscell(options.substrateExRxns)
        options.substrateExRxns = {options.substrateExRxns};
    end
    if ~isfield(options,'isAer')
        options.isAer = 1;
    end
    if ~isfield(options,'thKO')
        options.thKO = 0;
    end
    if ~isfield(options,'knockoutRxns')
        options.knockoutRxns = {};
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Set up model bounds


    % turn off glucose
    model = changeRxnBounds(model, 'EX_glc(e)', 0, 'l')

    % Knock out THs
    if options.thKO
        model = changeRxnBounds(model, thRxns, 0, 'b');
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% filename & waitbar info

    filename = sprintf('optSwapE_%s_%s.mat',...
                       options.filename, datestr(now,'yy-mm-dd_HH_MM_SS'));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Run swaps

    sortedTargetRxns = model.rxns(ismember(model.rxns, options.targetRxns));
    allSolns = [{'aer','thKO','subs','KOs','dh'} sortedTargetRxns'];
    count = 0;

    %%%% LOOP aerobic, anaerobic
    if options.isAer == 2
        aerFrom = 0; aerTo = 1;
    else
        aerFrom = options.isAer; aerTo = options.isAer;
    end

    lSize = aerTo-aerFrom+1;
    for l = aerFrom:aerTo

        % l is 0 for anaerobic, 1 for aerobic
        modelA = changeRxnBounds(model, 'EX_o2(e)',...
                                 -l*options.O2Uptake, 'l');


        %%%% LOOP exchange rxns
        mSize = length(options.substrateExRxns);
        for m = 1:length(options.substrateExRxns)

            substrateRxn = options.substrateExRxns{m};
            modelAS = changeRxnBounds(modelA, substrateRxn, ...
                                      -options.exUptakeMax, 'l');



            theseKnockoutRxns = options.knockoutRxns;
            modelAS = changeRxnBounds(modelAS, theseKnockoutRxns, 0, 'b');


            % % turn on FHL for H2 simulations
            % if (strcmp(targetRxn,'EX_h2(e)'))
            %     modelAS = changeRxnBounds(modelAS, 'FHL', 1000, 'u');
            % end

            %%%%%%%%%%
            % SINGLE Swaps
            %%%%%%%%%%

            %%%% LOOP dh rxns
            pSize = length(options.dhRxns)+1;
            for p = 0:length(options.dhRxns)

                if p < 1
                    newRxn = 'no swap';
                    modelASD = modelAS;
                else
                    % perform swap
                    [modelASD,newRxn] = ...
                        modelSwap(modelAS,options.dhRxns{p});
                end


                soln = optimizeCbModel(modelASD);
                display(['uptake ' ...
                         num2str(soln.x(ismember(modelASD.rxns,substrateRxn)))]);

                solnCell = {};
                solnFlux = soln.x(ismember(model.rxns,options.targetRxns));
                solnCell = num2cell(solnFlux)';
                allSolns = [allSolns; ...
                            {l,options.thKO,substrateRxn, theseKnockoutRxns,newRxn}, solnCell];


                count = count+1;
                display(sprintf('---- %g of %g', ...
                                count, lSize*mSize*(pSize + pSize^2)));
            end

            %%%%%%%%%%
            % DOUBLE Swaps
            %%%%%%%%%%

            display('double swaps')

            %%%% LOOP dh rxns
            pSize = length(options.dhRxns)+1;
            for p = 1:length(options.dhRxns)
                for q = 1:length(options.dhRxns)
                    modelASD = modelAS;

                    % perform swap
                    [modelASD,newRxn1] = ...
                        modelSwap(modelASD,options.dhRxns{p});
                    [modelASD,newRxn2] = ...
                        modelSwap(modelASD,options.dhRxns{q});
                    newRxnStr = [newRxn1 ' ' newRxn2];

                    soln = optimizeCbModel(modelASD);

                    solnCell = {};
                    solnFlux = soln.x(ismember(model.rxns,options.targetRxns));
                    solnCell = num2cell(solnFlux)';
                    allSolns = [allSolns; ...
                                {l,options.thKO,substrateRxn, theseKnockoutRxns,newRxnStr}, solnCell];
                    count = count+1;
                    display(sprintf('---- %g of %g', ...
                                    count, lSize*mSize*(pSize + pSize^2)));
                end
            end


        end


    end

    save(filename, 'allSolns')

end