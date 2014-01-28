function out = nonuniqueSwapsEcoli()

    % targets = {'EX_g3pg(e)'};
    targets = returnTargetRxns();
    dhs = dhRxnList('yield');

    substrates = {'EX_glc(e)', 'EX_xyl-D(e)', 'EX_glc(e)', 'EX_xyl-D(e)'};
    aerobicity = {'anaerobic', 'anaerobic', 'aerobic', 'aerobic'};
    out = cell(length(targets), length(substrates)*2);

    for i=1:length(substrates)
        [model, biomass] = setupModel('iJO', substrates(i), aerobicity(i));
        for k=1:length(targets)
            unique = {}; unique2 = {};
            max = 0;
            for l=1:length(dhs)
                temp = modelSwap(model, dhs(l), false);
                temp = changeRxnBounds(temp, biomass, 0.1, 'l');
                temp = changeObjective(temp, targets(k));
                soln = optimizeCbModel(temp);
                if (soln.f > 1.01*max)
                    max = soln.f;
                    unique = dhs(l);
                elseif (soln.f > 0.99*max)
                    unique{end+1} = dhs{l};
                end
            end
            % second swap
            for l=1:length(dhs)
                temp = modelSwap(model, unique{1}, false);
                temp = modelSwap(temp, dhs(l), false);
                temp = changeRxnBounds(temp, biomass, 0.1, 'l');
                temp = changeObjective(temp, targets(k));
                soln = optimizeCbModel(temp);
                if (soln.f > 1.01*max)
                    max = soln.f;
                    unique2 = dhs(l);
                elseif (soln.f > 0.99*max)
                    unique2{end+1} = dhs{l};
                end
            end
            out{k,i*2-1} = unique;
            out{k,i*2} = unique2;
        end
    end

    % save
    id = fopen('nonunique-ecoli-2.tsv', 'w');
    for i=1:size(out,1)
        for j=1:size(out, 2)
            for k=1:length(out{i,j})
                fprintf(id, '%s;', out{i,j}{k});
            end; 
            fprintf(id, '\t'); 
        end 
        fprintf(id, '\n'); 
    end
    fclose(id);
end

function rxns = returnTargetRxns() 
    rxns =  {'EX_cys-L(e)';
             'EX_hom-L(e)';
             'EX_thr-L(e)';
             'EX_15dap(e)';
             'EX_lys-L(e)';
             'EX_ile-L(e)';
             'EX_leu-L(e)';
             'EX_gly(e)';
             'EX_5mtr(e)';
             'EX_spmd(e)';
             'EX_pro-L(e)';
             'EX_cgly(e)';
             'EX_asp-L(e)';
             'EX_gthrd(e)';
             'EX_ser-L(e)';
             'EX_thym(e)';
             'EX_ptrc(e)';
             'EX_glyc(e)';
             'EX_pyr(e)';
             'EX_orn(e)';
             'EX_agm(e)';
             'EX_urea(e)';
             'EX_arg-L(e)';
             'EX_glyald(e)';
             'EX_lac-L(e)';
             'EX_asn-L(e)';
             'EX_thymd(e)';
             'EX_LalaDgluMdapDala(e)';
             'EX_acser(e)';
             'EX_ac(e)';
             'EX_LalaDgluMdap(e)';
             'EX_phe-L(e)';
             'EX_tyr-L(e)';
             'EX_ala-D(e)';
             'EX_alaala(e)';
             'EX_12ppd-S(e)';
             'EX_his-L(e)';
             'EX_g3pe(e)';
             'EX_g3pg(e)';
             'EX_enter(e)';
             'EX_feenter(e)';
             'EX_4abut(e)';
             'EX_ura(e)';
             'EX_ade(e)';
             'EX_pheme(e)'};
end