function productionOptSwapRKsingle

    count = '4';
    swaps = {'ACALD';};
    kos = {'EDA';'F6PA';'TPI';};
    kos_rk = {'F6PA';'FBA3';'G6PDH2r';'PFK';};
    substrate = 'EX_glc(e)';
    aerobicity = 'anaerobic';
    target = 'EX_ac(e)';

    filename = sprintf('c%s-%s-sub%s-%s',...
                       count,aerobicity,substrate,target);

    [model,biomass] = setupModel('iJO',substrate,aerobicity,'nothko');
    model_thko = setupModel('iJO',substrate,aerobicity,'thko');
    figure()
    hold all
    model_ko = changeRxnBounds(model, kos, 0, 'b')
    model_thko_ko = changeRxnBounds(model_thko, kos, 0, 'b')
    model_ko_swap = modelSwap(model_ko, swaps, false);
    model_thko_ko_swap = modelSwap(model_thko_ko, swaps, false);
    model_thko_rk = changeRxnBounds(model_thko,kos_rk,0,'b')
    model_rk = changeRxnBounds(model,kos_rk,0,'b')

    productionEnvelope(model,[],'k',target,biomass)
    productionEnvelope(model_thko,[],'--k',target,biomass)
    productionEnvelope(model_ko_swap,[],'r',target,biomass)
    productionEnvelope(model_thko_ko_swap,[],'--r',target,biomass)
    productionEnvelope(model_rk,[],'b',target,biomass)
    productionEnvelope(model_thko_rk,[],'--b',target,biomass)
    legend('wt','thko','design, no thko','design, thko',...
           'robustknock, nothko','robustknock, thko','Location','NorthEastOutside');
    title(filename,'Interpreter','None');
    set(gcf,'Color','White')
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf, 'Units', 'pixels');
    set(gcf, 'Position', [500 500 450 200]);
    % print('-dpng', filename);
    % close all


end
