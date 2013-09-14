function [dhRxns, matches] = locateDHs(model, mets)

    if nargin < 2
        nadh_met = 'nadh[c]';
        nad_met = 'nad[c]';
        nadph_met = 'nadph[c]';
        nadp_met = 'nadp[c]';
    else
        nadh_met = mets.nadh;
        nadph_met = mets.nadph;
        nad_met = mets.nad;
        nadp_met = mets.nadp;
    end
    theMet = nadh_met;

    nadhRxns = model.rxns(ismember(model.rxns,findRxnsFromMets(model,theMet)));
    nadhNames = model.rxnNames(ismember(model.rxns,nadhRxns));
    nadhSubsystems =  model.subSystems(ismember(model.rxns,nadhRxns));

    theMet = nadph_met;

    nadphRxns =  model.rxns(ismember(model.rxns,findRxnsFromMets(model,theMet)));
    nadphNames = model.rxnNames(ismember(model.rxns,nadphRxns));
    nadphSubsystems =  model.subSystems(ismember(model.rxns,nadphRxns));

    % Find NADH/NADPH reaction loops
    matches = {};

    % all NADH reactions
    S = full(model.S);
    rxnForms = S(:,find(ismember(model.rxns,nadhRxns)));
    %clear nadh/nad stoichiometry
    rxnForms(find(ismember(model.mets,nadh_met)),:) = zeros(1,length(nadhRxns));
    rxnForms(find(ismember(model.mets,nad_met)),:)  = zeros(1,length(nadhRxns));

    % all NADPH reactions
    pRxnForms = S(:,find(ismember(model.rxns,nadphRxns)));
    %clear nadh/nad stoichiometry
    pRxnForms(find(ismember(model.mets,nadph_met)),:) = zeros(1,length(nadphRxns));
    pRxnForms(find(ismember(model.mets,nadp_met)),:)  = zeros(1,length(nadphRxns));

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

    dhRxns = unique([nadhRxns; nadphRxns]);
end

