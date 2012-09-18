function [returnRxns,fluxes,soln] = almaasDistribution(model,options)

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

    dhRxns = locateDHs(model);



    % limit ammonia uptake rate to 100 mmol/g DW/h
    % model = changeRxnBounds(model, 'EX_nh4(e)', -100, 'l');

    % aerobic
    model = changeRxnBounds(model, 'EX_o2(e)', -20, 'l');

    % turn off glucose
    model = changeRxnBounds(model, 'EX_glc(e)', 0, 'l');

    % make possible loop reactions irreversible
    if ~isempty(length(options.possibleLoopRxns))
        model.rev(ismember(model.rxns,options.possibleLoopRxns)) = 0;
        model.lb(ismember(model.rxns,options.possibleLoopRxns)) = 0;
    end

    % h = waitbar(0,'almaasDistribution');



    fluxes = zeros(length(model.rxns),length(options.subs));
    for i = 1:length(options.subs)

        % set carbon substrate bounds
        modelTemp = changeRxnBounds(model, options.subs{i}, -8, 'l');

        if options.autRemLoops
            warning('this doesn''t work');
            if options.usePFBA
                error('cannot automatically remove loops with pFBA');
                return;
            end
            % use FBA
            soln = optimizeCbModel(modelTemp);
            remRxns = modelTemp.rxns(soln.x==-1000);
            modelTemp.rev(ismember(modelTemp.rxns,remRxns)) = 0;
            modelTemp.lb(ismember(modelTemp.rxns,remRxns)) = 0;
            
        end

        

        if options.usePFBA
            % use pFBA
            [~,~,modelIrrevFM] = pFBA(modelTemp,'skipclass',1);
            soln = optimizeCbModel(modelIrrevFM,'min');


            % modelIrrevFM = removeRxns(modelIrrevFM,'netFlux',false,false);
            % modelRev = convertToReversible(modelIrrevFM);
            % soln = optimizeCbModel(modelRev);
            % fluxesRev = soln.x;
            % % fluxesRev(ismember(modelRev.rxns,'netFlux')) = [];
            % fluxes(:,i) = fluxesRev;

            fluxesIrrev = soln.x;
            % remove pseudo reaction
            fluxesIrrev(ismember(modelIrrevFM.rxns,'netFlux')) = [];


            % convert irreversible fluxes to reversible ones
            [modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(modelTemp);
            selPosIrrevRxns = fluxesIrrev~=0;
            fluxesRev = zeros(length(modelTemp.rxns),1);
            fluxesRev(irrev2rev(selPosIrrevRxns)) = fluxesIrrev(selPosIrrevRxns);
            fluxes(:,i) = fluxesRev;

        else
            % use FBA
            soln = optimizeCbModel(modelTemp);
            fluxes(:,i) = soln.x;
        end

        if options.showPlot
            figure(i)
            % plot DHs
            % y = abs(soln.x(ismember(modelTemp.rxns,dhRxns)))./1000;
            % plot all reactions
            y = abs(soln.x)./1000;
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

        % waitbar(i/length(options.subs),h);
        display(sprintf('Optimization %g of %g', i, length(options.subs)));
    end


    % average fluxes over reactions
    avgFlux = mean(fluxes,2);
    % display([fluxes avgFlux])

    sel = abs(avgFlux)>options.threshFlux &...
          ismember(modelTemp.rxns,dhRxns);
    returnRxns = model.rxns(sel);
    fluxes = fluxes(sel);

    % save_to_base(1);

    % close(h)

    % display(size(returnRxns))
end