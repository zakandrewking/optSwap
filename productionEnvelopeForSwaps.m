function productionEnvelopeForSwaps(targetRxns,swaps,kos)
% productionEnvelopeForSwaps
%
% INPUTS
% swaps
% kos
% targetRxn
% keepWtDh
%
% Zachary King 9/12/12

    % cleaner = onCleanup(@() cleanup);
    % global run status
    status = 'starting';
    run = 'pipe test';
    
    % check inputs
    if nargin < 1
        targetRxns = {'EX_etoh(e)',
                      'EX_for(e)',
                      'EX_succ(e)',
                      'EX_ac(e)',
                      'EX_etoh(e)',
                      'EX_for(e)',
                      'EX_succ(e)',
                      'EX_ac(e)'}
          % targetRxn = 'EX_etoh(e)';

    end
    if nargin < 2
        swaps = {
            {}
            {}
            {}
            {}
            {'3OAR40';'3OAR80';'GLUDy';'LCARR';'ME1';}
            {'3OAR40';'EAR60x';'G6PDH2r';'GLUDy';'TRSARr';}
            {'G6PDH2r';'LCARR';'MDH';'MTHFD';'PDH';}
            {'GLUDy';'GND';'LCARR';'MDH';'PGCD';}
                };
    end
    if nargin < 3
        kos = {
            {'EDA';'F6PA';'FBA';'FLDR2';'TALA';}
            {'EDA';'F6PA';'PFK';'PGI';'TALA';}
            {'DHAPT';'EDA';'PFL';'PGI';'PYK';}
            {'EDA';'F6PA';'GHMT2r';'PGI';'TPI';}
            {}
            {}
            {}
            {}
              };
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
    [model, biomassRxn] = setupModel('iJO','EX_glc(e)',false,true);

    figure()
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
        
        printCoupledYield(modelT, targetRxns{i});
    end
    legendH = legend;
    set(legendH, 'Interpreter', 'none');
    ylabel('Production (mmol/gDW h)');
    xlabel('Growth rate (1/h)');
    hold off
    % title();
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