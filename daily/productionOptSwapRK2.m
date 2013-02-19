function productionOptSwapRK2

    f = fopen('/Users/zaking/Dropbox/lab/optSwap/data/optswap-robustknock-better.tsv', 'r');
    c = textscan(f,'%s%s%s%s%s%s%s','delimiter','\t');

    for i=1:length(c)
        labels(i) = c{i}(1);
    end
    script = c{1}(2:end);
    count = c{2}(2:end);
    target= c{3}(2:end);
    kos = parse(c{4}(2:end));
    swaps = parse(c{5}(2:end));
    aerobicity = c{6}(2:end);
    substrate = c{7}(2:end);

    [~,biomass] = setupModel('iJO','EX_glc(e)','anaerobic');

    for i=1:2:length(script)
        if strcmp(aerobicity{i},'anaerobic')
            continue
        end
        if ~strcmp(script{i},'OptSwap') || ~strcmp(script{i+1},'RobustKnock')
            error('bad order')
        end
        filename = sprintf('c%s-%s-sub%s-%s',...
                           count{i},aerobicity{i},substrate{i},target{i});

        model = setupModel('iJO',substrate{i},aerobicity{i},'nothko');
        model_thko = setupModel('iJO',substrate{i},aerobicity{i},'thko');
        % model = setupModelForTarget(model,target{i});
        % model_thko = setupModelForTarget(model_thko,target{i});
        
        figure()
        hold all
        model_ko = changeRxnBounds(model, kos{i}, 0, 'b');
        model_thko_ko = changeRxnBounds(model_thko, kos{i}, 0, 'b');
        model_ko_swap = modelSwap(model_ko, swaps{i}, false);
        model_thko_ko_swap = modelSwap(model_thko_ko, swaps{i}, false);
        model_thko_rk = changeRxnBounds(model_thko,kos{i+1},0,'b')
        model_rk = changeRxnBounds(model,kos{i+1},0,'b')

        productionEnvelope(model,[],'k',target{i},biomass)
        productionEnvelope(model_thko,[],'--k',target{i},biomass)
        productionEnvelope(model_ko_swap,[],'b',target{i},biomass)
        productionEnvelope(model_thko_ko_swap,[],'--b',target{i},biomass)
        productionEnvelope(model_rk,[],'r',target{i},biomass)
        productionEnvelope(model_thko_rk,[],'--r',target{i},biomass)
        legend('wt','thko','OptSwap, no thko','OptSwap, thko',...
               'RobustKnock, nothko','RobustKnock, thko','Location','NorthEastOutside');
        title(filename,'Interpreter','None');
        set(gcf,'Color','White')
        set(gcf, 'PaperPositionMode', 'auto');
        set(gcf, 'Units', 'pixels');
        set(gcf, 'Position', [500 500 450 200]);
        xlabel('growth rate (h^{-1})');
        ylabel('production (mmol gDW^{-1} h^{-1}');
        if strcmp(aerobicity{i},'aerobic')
            xlim([0 1.6]);
            set(gca,'XTick',0:0.2:2.0); 
        end
        print('-depsc2', filename)
        close all
    end
end

function theCell = parse(theCell)
    for i=1:length(theCell)
        if ~isempty(theCell{i})
            b = get_tokens(theCell{i}, ''';');
            b(strcmp('',b)) = [];
            theCell{i} = b;
        end
    end
end