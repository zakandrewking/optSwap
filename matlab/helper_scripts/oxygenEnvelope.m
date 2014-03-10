function oxygenEnvelope(model,targetRxn,options)
% Zachary King 1/31/2013

    if ~exist('options','var')
        options = struct();
    end
    if ~isfield(options, 'maxOxygenUptake'), options.maxOxygenUptake = 20; end
    if ~isfield(options, 'biomassRxn')
        options.biomassRxn = 'Ec_biomass_iJO1366_core_53p95M';
    end
    if ~isfield(options, 'nPts'), options.nPts = 20; end
    % if ~isfield(options, 'marker'), options.marker = []; end

    maxOxygenUptake = options.maxOxygenUptake;
    biomassRxn = options.biomassRxn;
    nPts = options.nPts;
    % marker = options.marker;

    % Create oxygen range vector
    oxygenValues = linspace(0,maxOxygenUptake,nPts);

    % set min biomass
    model = changeRxnBounds(model,biomass,0.1,'l');

    % Max/min for target production
    model = changeObjective(model,targetRxn);
    for i = 1:length(oxygenValues)
        model = changeRxnBounds(model,'EX_o2(e)',-oxygenValues(i),'l');
        sol = optimizeCbModel(model,'max');
        if (sol.stat > 0)
            targetUpperBound(i) = sol.f;
        else
            targetUpperBound(i) = NaN;
        end
        sol = optimizeCbModel(model,'min');
        if (sol.stat > 0)
            targetLowerBound(i) = sol.f;
        else
            targetLowerBound(i) = NaN;
        end
    end

    % Plot results
    lineHandle=plot([oxygenValues fliplr(oxygenValues)],[targetUpperBound fliplr(targetLowerBound)],'LineWidth',2);
    axis tight;
    ylabel([strrep(targetRxn,'_','-') ' (mmol/gDW h)']);
    xlabel('oxygen uptake rate (mmol/gDW-h)');
    set(gcf,'Color','White'); 
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf, 'Units', 'pixels');
    set(gcf, 'Position', [500 500 400 280]);
end