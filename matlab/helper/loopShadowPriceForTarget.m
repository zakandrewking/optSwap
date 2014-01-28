function loopShadowPriceForTarget(targets)

    file = 'shadow_prices.txt';
    fileId = fopen(file, 'a');
    % fileId = [];
    [model, biomass] = setupModel('iJO', 'EX_glc(e)', 'anaerobic', 'thko');
    
    for i = 1:length(targets)
        [mets, prices] = shadowPriceForTarget(model, biomass, targets{i});
        if i==1
            fprintf(fileId, '\t');
            for j = 1:length(mets)
                fprintf(fileId, '%s\t', mets{j});
            end
            fprintf(fileId, '\n');
        end
        fprintf(fileId, '%s\t', targets{i});
        for j = 1:length(prices)
            fprintf(fileId, '%f\t', prices(j));
        end
        fprintf(fileId, '\n');
    end
    fclose(fileId);
end