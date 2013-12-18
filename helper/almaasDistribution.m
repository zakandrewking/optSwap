function [returnRxns, fluxes, f] = almaasDistribution(model,options)

% ALMAASDISTRIBUTION
%
% Find distribution of fluxes after FBA of model for biomass
% on glutamate and succinate, as described in [1].
% (aerobic conditions)
%
% INPUTS
% model - cobra model
% options
%   subs - substrates to consider
%   threshFlux - threshold
%   showPlot
%   usePFBA - true: usePFBA, false: FBA (default)
%   possibleLoopRxns -
%   autRemLoops - true: automatically remove possible loops, false (default)
%   dhCount
%
%
% REF
%  1. Almaas, E., Kovács, B., Vicsek, T., Oltvai, Z. N.,
%     & Barabási, a-L. (2004). Global organization of
%     metabolic fluxes in the bacterium Escherichia coli.
%     Nature, 427(6977), 839–43. doi:10.1038/nature02289
%
% Zachary King, 7/1/2012
%

    if nargin < 1, error('Not enough arguments'); end
    if ~exist('options','var')
        options = struct();
    end
    if ~isfield(options,'subs'), options.subs = {'EX_glc(e)'}; end
    if ~iscell(options.subs), options.subs = {options.subs}; end

    if ~isfield(options,'usePFBA'), options.usePFBA = false; end
    if ~isfield(options,'threshFlux'), options.threshFlux = 0.2; end
    if ~isfield(options,'showPlot'), options.showPlot = true; end
    if ~isfield(options,'possibleLoopRxns'), options.possibleLoopRxns = {}; end
    if ~isfield(options,'autRemLoops'), options.autRemLoops = false; end
    if ~isfield(options,'glucoseExchange'), options.glucoseExchange = 'EX_glc(e)'; end
    if ~isfield(options,'exchangeMag'), options.exchangeMag = 10; end

    if isfield(options, 'metStruct')
        dhRxns = locateDHs(model, options.metStruct);
    else
        dhRxns = locateDHs(model);
    end
    fprintf('found %d oxidoreductase reactions\n', length(dhRxns));

    % turn off glucose
    model = changeRxnBounds(model, options.glucoseExchange, 0, 'l');

    % make possible loop reactions irreversible
    if ~isempty(length(options.possibleLoopRxns))
        model.rev(ismember(model.rxns,options.possibleLoopRxns)) = 0;
        model.lb(ismember(model.rxns,options.possibleLoopRxns)) = 0;
    end

    fluxes = zeros(length(model.rxns),length(options.subs));
    for i = 1:length(options.subs)
        % waitbar(i/length(options.subs),h);

        % set carbon substrate bounds
        modelTemp = changeRxnBounds(model, options.subs{i}, -options.exchangeMag, 'l');

        % use FBA
        soln_fba = optimizeCbModel(modelTemp);
        f = soln_fba.f;
        if f == 0
            display('no growth')
            fluxes(:,i) = zeros(size(fluxes(:,1)));
        else
            % only run pFBA if the cell is viable
            if options.usePFBA
                % use pFBA
                [~,~,modelIrrevFM] = pFBA(modelTemp,'geneoption',0,'tol',1e-7,'skipclass',1);
                soln_pfba = optimizeCbModel(modelIrrevFM);

                fluxesIrrev = soln_pfba.x;
                % remove pseudo reaction
                fluxesIrrev(ismember(modelIrrevFM.rxns,'netFlux')) = [];

                % convert irreversible fluxes to reversible ones
                [modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(modelTemp);
                selPosIrrevRxns = fluxesIrrev~=0;
                fluxesRev = zeros(length(modelTemp.rxns),1);
                fluxesRev(irrev2rev(selPosIrrevRxns)) = fluxesIrrev(selPosIrrevRxns);
                fluxes(:,i) = fluxesRev;
            else
                fluxes(:,i) = soln_fba.x;
            end
        end

        if options.showPlot
            figure(i)
            % plot DHs
            % y = abs(soln.x(ismember(modelTemp.rxns,dhRxns)))./1000;
            % plot all reactions
            y = abs(fluxes(:,i))./1000;
            x = [];
            for j=-5:.2:0
                x(end+1) = 10^(j);
            end
            n = histc(y,x)./length(y);
            plot(x,n,'o');
            set(gca,'xscale','log')
            set(gca,'yscale','log')
            title(options.subs{i},'Interpreter', 'none')
        end
    end

    % average fluxes over reactions
    avgFlux = mean(fluxes,2);
    % just select fluxes for dh reactions
    avgFlux = avgFlux(ismember(model.rxns, dhRxns));
    % find the largest absolute fluxes
    avgFlux = abs(avgFlux);
    % keep track of indices
    avgFlux(:,2) = (1:length(avgFlux))';
    % sort descending by column 1
    avgFlux = sortrows(avgFlux,-1);
    count = min(length(dhRxns), options.dhCount);
    % initialize return reactions
    returnRxns = cell(count,1);
    for i=1:length(returnRxns)
        returnRxns(i) = dhRxns(avgFlux(i,2));
    end
    fluxes = avgFlux(1:count,1);

    % fprintf('flux through (R,R)-butanediol dehydrogenase: %.2f\n', ...
    %         soln_fba.x(ismember(modelTemp.rxnNames, ['(R,R)-butanediol ' ...
    %                     'dehydrogenase'])));
    % model.rxnNames(ismember(model.rxns,returnRxns{1}))
    % fluxes(1)
    % fprintf('flux through glyceraldehyde-3-phosphate dehydrogenase: %.2f\n', ...
    %         soln_fba.x(ismember(modelTemp.rxnNames, 'glyceraldehyde-3-phosphate dehydrogenase')));    
end