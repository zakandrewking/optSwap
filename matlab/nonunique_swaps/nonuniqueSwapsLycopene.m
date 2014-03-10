function out = nonuniqueSwapsLycopene()

    targets = {'EX_lyco(e)'};
    dhs = dhRxnList('yield');

    substrates = {'EX_glc(e)', 'EX_xyl_D(e)', 'EX_glyc(e)', ...
                  'EX_glc(e)', 'EX_xyl_D(e)', 'EX_glyc(e)'};
    aer = {'anaerobic','anaerobic','anaerobic','aerobic','aerobic','aerobic'};
    out = cell(length(targets), length(substrates)*2);
    maxes = zeros(size(out));

    for i=1:length(substrates)
        [model, biomass] = setupModel('iJO-h', substrates(i), aer(i));
        model = makeLycopene(model);
        for k=1:length(targets);
            unique = {}; unique2 = {};
            max = 0; max2 = 0;
            for l=1:length(dhs)
                temp = modelSwap(model, dhs(l), false);
                if (i==3)
                    s_wt = optimizeCbModel(model);
                    min_biomass = 0.1*s_wt.f;
                else
                    min_biomass = 0.1;
                end
                temp = changeRxnBounds(temp, biomass, min_biomass, 'l');
                temp = changeObjective(temp, targets(k));
                soln = optimizeCbModel(temp);
                if (soln.f > 1.01*max)
                    max = soln.f;
                    unique = dhs(l);
                elseif (soln.f > 0.99*max)
                    unique{end+1} = dhs{l};
                end
            end
          
            if length(unique) > 0
                % second swap
                for l=1:length(dhs)
                    temp = modelSwap(model, unique{1}, false);
                    temp = modelSwap(temp, dhs(l), false);
                    temp = changeRxnBounds(temp, biomass, min_biomass, 'l');
                    temp = changeObjective(temp, targets(k));
                    soln = optimizeCbModel(temp);
                    if (soln.f > 1.01*max2)
                        max2 = soln.f;
                        unique2 = dhs(l);
                    elseif (soln.f > 0.99*max2)
                        unique2{end+1} = dhs{l};
                    end
                end
            end
            out{k,i*2-1} = unique;
            out{k,i*2} = unique2;
            maxes(k,i*2-1) = max;
            maxes(k,i*2) = max2;
        end
    end

    % save
    id = fopen('nonunique-lycopeen.tsv', 'w');
    for i=1:size(out,1)
        for j=1:size(out, 2)
            for k=1:length(out{i,j})
                fprintf(id, '%s;', out{i,j}{k});
            end; 
            fprintf(id, '\t'); 
        end 
        fprintf(id, '\n'); 
    end
    fprintf(id, '\n'); 
    for i=1:size(maxes,1)
        for j=1:size(maxes, 2)
            for k=1:length(maxes(i,j))
                fprintf(id, '%f', maxes(i,j));
            end; 
            fprintf(id, '\t'); 
        end 
        fprintf(id, '\n'); 
    end
    fclose(id);
end
