function model = turnOnSubsystem(model, m, subSystems)
    these_rxns_irrev = model.rxns(ismember(model.subSystems,subSystems{m}) ...
                                  & model.rev==0);
    these_rxns_rev = model.rxns(ismember(model.subSystems,subSystems{m}) ...
                                & model.rev==1);
    if ~isempty(these_rxns_irrev)
        model = changeRxnBounds(model,these_rxns_irrev,1000, 'u');
    end
    if ~isempty(these_rxns_rev)
        model = changeRxnBounds(model,these_rxns_rev,  1000,'u');
        model = changeRxnBounds(model,these_rxns_rev, -1000, 'l');
    end
end