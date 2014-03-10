function productionEnvelopeForSwaps(modelName,targetRxns,swaps,kos,aerobicString,substrate)
% productionEnvelopeForSwaps
%
% INPUTS
% targetRxn
% swaps
% kos
% aerobicString
% substrate
%
% Zachary King 9/12/12

% cleaner = onCleanup(@() cleanup);
% global run status
    status = 'starting';
    run = 'pipe test';

    % check inputs 
    if nargin < 1
        targetRxns = 'EX_ac(e)';
    end
    if nargin < 2
        swaps = {
            {}
            {'GLUDy';'MDH';'PGCD';}
            {}
            {'GLUDy';'MDH';'PGCD';}
            {'GLUDy';'GND';'LCARR';'MDH';'PGCD'}
                };
    end
    if nargin < 3
        kos = {
            {}
            {}
            {'DRPA';'F6PA';'TPI';}
            {'EDA';'F6PA';'TPI';}
            {'EDA';'F6PA';'PGM';'PPKr';'TPI';}
              };
    end
    if nargin < 4
        aerobicString = 'anaerobic';
    end  
    if nargin < 5
        substrate = 'EX_glc(e)';
    end
      
    if isempty(swaps)
        swaps = cell(size(kos));
    elseif isempty(kos)
        kos = cell(size(swaps));
    end
    if ~iscell(targetRxns)
        aTargetRxn = targetRxns;
        targetRxns = {};
        for i=1:max(length(swaps),length(kos))
            targetRxns = [targetRxns; {aTargetRxn}];
        end
    end

    % setup model
    keepWtDh = false;
    global biomassRxn
    [model, biomassRxn] = setupModel(modelName,substrate,aerobicString,'nothko');

    color = {'k', 'Blue', 'Red', 'Green'};
    lineStyle = {'-', '--', '-', '--'};
    
    %figure()
    hold on
    for i=1:length(targetRxns)
        status = sprintf('run%d', i);
        modelT = model;
        modelT = setupModelForTarget(modelT, targetRxns{i});
        if ~isempty(swaps{i})
            modelT = modelSwap(modelT, swaps{i}, keepWtDh);
        end
        if ~isempty(kos{i})
            modelT = changeRxnBounds(modelT, kos{i}, 0, 'b')
        end
        [~,~,lineHandle] = ...
            productionEnvelope(modelT, [], 'k', targetRxns{i}, biomassRxn);
        koStr='';
        for j=1:length(kos{i})
            koStr = [koStr ' ' kos{i}{j}];
        end
        swapStr='';
        for j=1:length(swaps{i})
            swapStr = [swapStr ' ' swaps{i}{j}];
        end
        set(lineHandle,'DisplayName',[koStr '; ' swapStr '--' ...
                            targetRxns{i}]);
        if i < length(color)
            set(lineHandle, 'Color', color{i});
            set(lineHandle, 'LineStyle', lineStyle{i});
            set(lineHandle, 'LineWidth', 1.5);
        end

        printCoupledYield(modelT, targetRxns{i});
    end
    hYLabel = ylabel('Production (mmol/gDW h)');
    hXLabel = xlabel('Growth rate (1/h)');
    
    hLegend = legend('show');
    set(hLegend, 'Interpreter', 'none');
    set([hLegend gca hXLabel, hYLabel], 'FontSize', 10);
    set([gca hXLabel, hYLabel], 'FontName', 'Helvetica');
    set(gcf, 'Color', 'White')
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf, 'Units', 'pixels');
    set(gcf, 'Position', [500 500 200 140]);
    set(gca, ...
        'Box'         , 'on'     , ...
        'TickDir'     , 'in'     , ...
        'TickLength'  , [.015 .015] , ...
        'XMinorTick'  , 'off'      , ...
        'YMinorTick'  , 'off'      , ...
        'YTick'       , 0:10:80, ...
        'XTick', 0:0.1:2.0, ... %'XTick', 0:0.2:2.0, ...
        'LineWidth'   , 1         );
    
    %set(hTitle, 'FontSize', 12, 'FontWeight', 'bold');
    legend('off')
    hold off
    print -depsc2 finalPlot.eps
    status = 'finished';
end


function printCoupledYield(model, targetRxn)
    global biomassRxn
    model.c = zeros(size(model.c));
    model.c(ismember(model.rxns, biomassRxn)) = 1;
    soln = optimizeCbModel(model);
    if ~isempty(soln.x)
        display(sprintf('%.4f,', soln.x(ismember(model.rxns, targetRxn))));
    else
        display('no sol,');
    end
end