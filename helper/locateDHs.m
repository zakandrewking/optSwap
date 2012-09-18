function [dhRxns, matches] = locateDHs(model)

    theMet = 'nadh[c]';

    nadhRxns = model.rxns(ismember(model.rxns,findRxnsFromMets(model,theMet)));
    nadhNames = model.rxnNames(ismember(model.rxns,nadhRxns));
    nadhSubsystems = ...
        model.subSystems(ismember(model.rxns,nadhRxns));
    % xlswrite('iJO-DHs.xls',rxnNames);

    theMet = 'nadph[c]';

    nadphRxns =  model.rxns(ismember(model.rxns,findRxnsFromMets(model,theMet)));
    nadphNames = model.rxnNames(ismember(model.rxns,nadphRxns));
    nadphSubsystems = ...
        model.subSystems(ismember(model.rxns,nadphRxns));
    % xlswrite('iJO-DHs.xls',rxnNames);

    % Find NADH/NADPH reaction loops
    matches = {};

    % all NADH reactions
    S = full(model.S);
    rxnForms = S(:,find(ismember(model.rxns,nadhRxns)));
    %clear nadh/nad stoichiometry
    rxnForms(find(ismember(model.mets,'nadh[c]')),:) = zeros(1,length(nadhRxns));
    rxnForms(find(ismember(model.mets,'nad[c]')),:)  = zeros(1,length(nadhRxns));

    % all NADPH reactions
    pRxnForms = S(:,find(ismember(model.rxns,nadphRxns)));
    %clear nadh/nad stoichiometry
    pRxnForms(find(ismember(model.mets,'nadph[c]')),:) = zeros(1,length(nadphRxns));
    pRxnForms(find(ismember(model.mets,'nadp[c]')),:)  = zeros(1,length(nadphRxns));

    for i = 1:length(nadhRxns)
        for j = 1:length(nadphRxns)
            if (isequal(rxnForms(:,i),pRxnForms(:,j)))
                matches = [matches; [nadhRxns(i), nadphRxns(j)]]; 
                if (strcmp(nadhRxns{i},'3HCINNMH'))
                    display('stop')
                end
                break;
            end
        end
    end

    dhRxns = [nadhRxns; nadphRxns];

    save_to_base(1);

end

