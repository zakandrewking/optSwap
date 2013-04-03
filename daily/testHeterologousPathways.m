exch_reactions = {'EX_13ppd(e)';
                  'EX_14btd(e)';
                  'EX_btd_RR(e)';
                  'EX_btd_SS(e)';
                  'EX_btd_meso(e)';
                  'EX_bhb(e)';
                  'EX_bhb_S(e)';
                  'EX_3hpp(e)';
                  'EX_R_3hpt(e)';
                  'EX_S_3hpt(e)';
                  'EX_styr(e)';
                  'EX_phe_L(e)';
                  'EX_phstyr(e)'
                 };

subsystems = {' 1,3-Propanediol production';
              ' 1,4-Butanediol production'; 
              ' 2,3-Butanediol production';
              ' 2,3-Butanediol production';
              ' 2,3-Butanediol production';
              ' 3-Hydroxybutyrate production';
              ' 3-Hydroxybutyrate production';
              ' 3-Hydroxypropanoate production'; 
              ' 3-Hydroxyvalerate production'; 
              ' 3-Hydroxyvalerate production'; 
              'Styrene production';
              'Styrene production';
              'p-Hydroxystyrene production'};

intermediate = {'Production pathway intermediate synthesis';};


substrates = {'EX_glc(e)', 'EX_xyl_D(e)'};
aer = {'anaerobic','aerobic'};
thko = {'thko', 'nothko'}; 
for i=1:length(substrates)
    for j=1:length(aer)
        for k=1:length(thko) 
            model = setupModel('iJO-h',substrates{i},aer{j},thko{k});
            for m = 1:length(exch_reactions)
                model_t = model;
                these_rxns_irrev = model.rxns(ismember(model.subSystems,subsystems{m}) ...
                                              & model.rev==0);
                these_rxns_rev = model.rxns(ismember(model.subSystems,subsystems{m}) ...
                                              & model.rev==1);
                if ~isempty(these_rxns_irrev)
                    model_t = changeRxnBounds(model_t,these_rxns_irrev,1000, 'u');
                end
                if ~isempty(these_rxns_rev)
                    model_t = changeRxnBounds(model_t,these_rxns_rev,  1000,'u');
                    model_t = changeRxnBounds(model_t,these_rxns_rev, -1000, 'l');
                end
                model_t = changeObjective(model_t, exch_reactions{m});
                soln = optimizeCbModel(model_t);
                fprintf('%s\t%s\t%s\t%s\t%f\n', substrates{i}, aer{j}, thko{k}, ...
                        exch_reactions{m}, soln.f);
            end 
        end
    end
end