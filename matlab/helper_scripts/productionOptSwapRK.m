function productionOptSwapRK

    f = fopen('/Users/zaking/Dropbox/lab/optSwap/data/optswap-robustknock-runs.tsv', 'r');
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

    if false
        filename = sprintf('c%s-%s-sub%s-%s-%s',...
                           count{i},aerobicity{i},substrate{i},target{i},script{i});

        model = setupModel('iJO',substrate{i},aerobicity{i},'nothko');
        model_thko = setupModel('iJO',substrate{i},aerobicity{i},'thko');

        figure()
        hold all
        model_ko = changeRxnBounds(model, kos{i}, 0, 'b');
        model_thko_ko = changeRxnBounds(model_thko, kos{i}, 0, 'b');
        model_ko_swap = modelSwap(model_ko, swaps{i}, false);
        model_thko_ko_swap = modelSwap(model_thko_ko, swaps{i}, false);
        % model_thko_rk = changeRxnBounds(model_thko,kos_rk,0,'b')
        % model_rk = changeRxnBounds(model,kos_rk,0,'b')

        productionEnvelope(model,[],'k',target{i},biomass)
        productionEnvelope(model_thko,[],'--k',target{i},biomass)
        productionEnvelope(model_ko_swap,[],'r',target{i},biomass)
        productionEnvelope(model_thko_ko_swap,[],'--r',target{i},biomass)
        % productionEnvelope(model_rk,[],'b',target,biomass)
        % productionEnvelope(model_thko_rk,[],'--b',target,biomass)
        legend('wt','thko','design, no thko','design, thko','Location','NorthEastOutside') %,...
               % 'robustknock, nothko','robustknock, thko');
        title(filename,'Interpreter','None');
        set(gcf,'Color','White')
        set(gcf, 'PaperPositionMode', 'auto');
        set(gcf, 'Units', 'pixels');
        set(gcf, 'Position', [500 500 450 200]);
        print('-dpng', filename);
        close all
    end
    
    for i=1:length(script)
        model_ko = changeRxnBounds(model, kos{i}, 0, 'b');
        model_ko_swap = modelSwap(model_ko, swaps{i}, false);
        soln = optimizeCbModel(model);
        if abs(soln.x(ismember(model.rxns,this_swap{j}))) < 0.01
            display('this swap is not necessary')
        end
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