function runAlmaasDistributionYeast

    count = 200;
    % glc, D-xyl, gylc, L-arab
    % subs = {'r_1714', 'r_1718', 'r_1808', 'r_1878', ...
    %         'r_1714', 'r_1718', 'r_1808', 'r_1878'};
    subs = {'EX_glc(e)', 'EX_xyl-D(e)', ...
            'EX_glc(e)', 'EX_xyl-D(e)'};
    aerString = {'aerobic', 'aerobic', ...
                 'anaerobic', 'anaerobic'};
    s = length(subs);
    out = cell(s,1);
    f_out = zeros(s,1);
    for i=1:s
        fprintf('Case %d of %d\n', i, length(subs));
        model = setupModel('iMM904',subs{i},aerString{i});
        model = changeRxnBounds(model, 'ALCD2ir', 0, 'b');
        options.subs = subs{i};
        options.usePFBA = true;
        options.dhCount = count;
        options.showPlot = false;
        % options.glucoseExchange = 'r_1714';
        options.glucoseExchange = 'EX_glc(e)';
        % options.metStruct = struct('nadh', 's_1203', 'nadph', 's_1212', ...
        %                            'nad', 's_1198', 'nadp', 's_1207');
        [returnRxns, fluxes, f] = almaasDistribution(model,options);
        out{i} = [returnRxns, num2cell(fluxes)];
        f_out(i) = f;
    end
    fileId = fopen('almaas_iMM904_all_dh_output_pFBA_ALCD2ir-ko.tsv', 'w');
    fprintf(fileId, 'glc\t\tD-xyl\t\tglc\t\tD-xyl\n');
    

    for j=0:count*2
        for i=1:s
            if (j==0)
                fprintf(fileId, 'growth rate\t%.3f\t', f_out(i));
            elseif (j <= size(out{i},1))
                r_name_cell = model.rxnNames(ismember(model.rxns, out{i}{j,1}));
                fprintf(fileId, '%s\t%f\t', r_name_cell{1}, out{i}{j,2});
            end
        end 
        fprintf(fileId, '\n');
    end
end