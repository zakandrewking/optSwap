function twoModels(model1,model2)
    
    size(model1.S)
    size(model2.S)
    
    find(model1.c)
    find(model2.c)
    
    for i = 1:length(model1.rxns)
        rxn = model1.rxns{i};
        if ~find(ismember(model2.rxns, rxn))
            fprintf('reaction %s missing from model 2\n\n')
        end
        lb1 = model1.lb(ismember(model1.rxns,rxn));
        lb2 = model2.lb(ismember(model2.rxns,rxn));
        ub1 = model1.ub(ismember(model1.rxns,rxn));
        ub2 = model2.ub(ismember(model2.rxns,rxn)); 
        rev1 = model1.rev(ismember(model1.rxns,rxn));
        rev2 = model2.rev(ismember(model2.rxns,rxn));
        if lb1 ~= lb2
            fprintf('%s lower bounds: %d / %d\n', rxn, lb1, lb2);
        end
        if ub1 ~= ub2
            fprintf('%s upper bounds: %d / %d\n', rxn, ub1, ub2);
        end
        if rev1 ~= rev2
            fprintf('%s reversibilities: %d / %d\n', rxn, rev1, rev2);
        end
    end
    
   
end