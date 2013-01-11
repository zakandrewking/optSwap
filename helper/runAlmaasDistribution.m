function runAlmaasDistribution

    subs = {'EX_glc(e)';'EX_glc(e)';'EX_xyl-D(e)';'EX_xyl-D(e)'};
    aerString = {'anaerobic','aerobic','anaerobic','aerobic'};
    for i=1:4
        model = setupModel('iJO','EX_glc(e)',aerString{i},'noTHKO');
        options.subs = subs{i};
        options.possibleLoopRxns = {'TRSARr'};
        % options.autRemLoops = true;
        % options.usePFBA = true;
        options.dhCount = 31;
        options.showPlot = false;
        [returnRxns,fluxes] = almaasDistribution(model,options);
        out{i} = [returnRxns, num2cell(fluxes)];
        % save_to_base(1);
        % S = tw.updateStatus('this gremlin is awake');
    end
    fileId = fopen('almaas_output.csv', 'w')
    label = {'glc anaerobic\n','glc aerobic\n','xyl anaerobic\n','xyl aerobic\n'};
    for j=1:4
        fprintf(fileId,label{j});
        for i=1:size(out{j},1)
            fprintf(fileId, '%s,%f\n',out{j}{i,1}, out{j}{i,2});
        end 
    end
end