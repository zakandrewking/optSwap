function soln = optSwapB(model,thRxns,parameterCell,options)

% OPTSWAPB
%
% Swaps electron carrier specificity for dehydrogenase reactions. Runs FBA
% on TH knockout model and unrestricted TH model. Accepts a cell with all
% parameter options for greater control.
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
% parameterCell - includes target reactions, th knockout, substrates,
%       knockouts, and aerobicity in one matrix. Form:
%       { isAer, thKO, subs, target, knockouts }
%
% OPTIONAL
% options.
%   dhRxns - List of dehydrogenase reactions to swap.
%   exUptakeMax - default 20
%   O2Uptake - default 20
%   filename
%
% OUTPUT
% soln
%
% Zachary King 7/6/2012


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check options & set defaults

    if nargin < 3, error('not enough arguments'); end

    if ~exist('options','var')
        options = struct();
    end
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

    %       { isAer, thKO, substrate, target, knockouts }
    parameterStruct.isAer     = parameterCell(:,1);
    parameterStruct.thKO      = parameterCell(:,2);
    parameterStruct.substrate = parameterCell(:,3);
    parameterStruct.target    = parameterCell(:,4);
    parameterStruct.knockouts = parameterCell(:,5);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Set up model bounds

    % turn off glucose
    model = changeRxnBounds(model, 'EX_glc(e)', 0, 'l');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% filename & waitbar info

    filename = sprintf('optSwapB_%s_%s.mat',...
                       options.filename, datestr(now,'yy-mm-dd_HH_MM_SS'));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Run swaps

    soln = {'aer','thKO','subs','target','KOs','dh','YPS','SSP','SOC'};
    count = 0;
    ticID = tic;
    
    iSize = size(parameterCell,1);
    for i = 1:iSize

        
        modelT = model;

        % set aerobicity
        isAer = parameterStruct.isAer{i}
        if isAer
            modelT = changeRxnBounds(modelT, 'EX_o2(e)',...
                                     -options.O2Uptake, 'l');
        else
            modelT = changeRxnBounds(modelT, 'EX_o2(e)', 0, 'l');
        end

        % Knock out THs
        thKO =  parameterStruct.thKO{i};
        if parameterStruct.thKO{i}
            modelT = changeRxnBounds(modelT, thRxns, 0, 'b');
        end

        % set substrate rxn
        substrateRxns = parameterStruct.substrate{i}
        modelT = changeRxnBounds(modelT, substrateRxns, ...
                                 -options.exUptakeMax, 'l');

        % set target rxn
        targetRxn = parameterStruct.target{i};

        % turn on FHL for H2 simulations
        if (strcmp(targetRxn,'EX_h2(e)'))
            modelT = changeRxnBounds(modelT, 'FHL', 1000, 'u');
        end

        % set knockouts
        knockoutRxns = parameterStruct.knockouts{i};
        modelT = changeRxnBounds(modelT, knockoutRxns, 0, 'b');

        pSize = length(options.dhRxns);
        for p = 0:pSize

            if p < 1
                newRxn = 'no swap';
                modelTD = modelT;
            else
                % perform swap
                [modelTD,newRxn] = ...
                    modelSwap(modelT,options.dhRxns{p});
            end

            % show production envelopes
            biomassRxn = model.rxns(model.c==1);
            [biomassValues,targetValues,~] = ...
                productionEnvelope(modelTD,knockoutRxns, 'k',...
                                   targetRxn,biomassRxn);
            figure(1)
            plot(biomassValues,targetValues);
            title(sprintf('%g %s %s %s',...
                          isAer,substrateRxns,targetRxn,newRxn));
            [~,ind] = max(biomassValues);
            YPS = targetValues(end,1)/options.exUptakeMax;
            SSP = YPS * biomassValues(end);
            lastZeroInd = find(targetValues(:,1),1,'first')-1;
            if lastZeroInd < 1
                SOC = [];
            else
                slope = ...
                    (targetValues(end,1) - targetValues(lastZeroInd,1))...
                    ./ ...
                    (biomassValues(end,1) - biomassValues(lastZeroInd,1));
                SOC = YPS.^2 ./ slope;
            end

            %display([num2str(i) ' ' thisNewRxn ' ' num2str(targetValues(ind))]);

            % opt = optimizeCbModel(modelASD);
            % display(['uptake ' num2str(opt.x(ismember(modelASD.rxns,substrateRxn)))]);
            % dhFlux = soln.x(ismember(thisModel.rxns,thisNewRxn));
            % display(['Flux through ' thisNewRxn ': ' num2str(dhFlux)]);

            if iscell(substrateRxns), substrateString = makeString(substrateRxns);
            else, substrateString = substrateRxns; end
            if iscell(knockoutRxns), knockoutString = makeString(knockoutRxns);
            else, knockoutString = knockoutRxns; end
            
            soln(end+1,:) = {isAer,thKO,substrateString,...
                             targetRxn,knockoutString,newRxn,...
                             YPS,SSP,SOC};

            save(filename, 'soln')
            % clc
            % TODO count if off
            count = count+1;
            display(sprintf('---- %g of %g', ...
                            count, iSize*(pSize+1)));

            if count<10
                t = toc(ticID);
                estTime = (t/count*1)*iSize*(pSize+1)/60;
                display(sprintf('time estimate: %g min', estTime));
            end
        end

    end
end

function aStr = makeString(cell)
    aStr = '';
    for i=1:size(cell,1)
        for j=1:size(cell,2)
            aStr = [aStr ' ' cell{i,j}];
        end
    end
end