function phenotypeSurface(model,targetRxn,options)

    if ~exist('options','var')
        options = struct();
    end
    if ~isfield(options, 'maxOxygenUptake'), options.maxOxygenUptake = 20; end
    if ~isfield(options, 'biomassRxn')
        options.biomassRxn = 'Ec_biomass_iJO1366_core_53p95M';
    end
    if ~isfield(options, 'nPts'), options.nPts = 20; end
    if ~isfield(options, 'faceColor'), options.faceColor = 'g'; end
    
    maxOxygenUptake = options.maxOxygenUptake;
    biomassRxn = options.biomassRxn;
    nPts = options.nPts;
    faceColor = options.faceColor;

    % Run FBA to get upper bound for biomass
    model = changeObjective(model,biomassRxn);
    model = changeRxnBounds(model,'EX_o2(e)',-maxOxygenUptake,'l');
    solMax = optimizeCbModel(model,'max');
    solMin = optimizeCbModel(model,'min');

    % Create biomass range vector
    biomassValues = linspace(solMin.f,solMax.f,nPts);

    % Create oxygen range vector
    oxygenValues = linspace(0,maxOxygenUptake,nPts);

    % Max/min for target production
    model = changeObjective(model,targetRxn);
    for i=1:length(oxygenValues)
        for j=1:length(biomassValues)
            % fprintf('o2 %f   ', -oxygenValues(i));
            % fprintf('biomass %f\n', biomassValues(j));
            model = changeRxnBounds(model,'EX_o2(e)',-oxygenValues(i),'b');
            model = changeRxnBounds(model,biomassRxn,biomassValues(j),'b');
            sol = optimizeCbModel(model);
            if (sol.stat > 0)
                yield(j,i) = sol.f;
            else
                yield(j,i) = 0;
            end
        end
    end

    % Plot results
    lineHandle=mesh(oxygenValues,biomassValues,yield,...
                    'FaceAlpha',0.7,'EdgeAlpha',0.7,'EdgeColor','k',...
                    'FaceColor',faceColor);
    zlabel([strrep(targetRxn,'_','-') ' (mmol gDW^{1} h^{-1})']);
    xlabel('oxygen uptake rate (mmol gDW^{-1} h^{-1})');
    ylabel('growth rate (h^{-1})');
    set(gcf,'Color','White'); 
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf, 'Units', 'pixels');
    set(gcf, 'Position', [500 500 700 500]);
    xlim([0 maxOxygenUptake]);
    ylim([solMin.f solMax.f]);
    hold on
end