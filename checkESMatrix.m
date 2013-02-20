function [isValid, rxnsFound] = checkESMatrix(model)
% check whether E*S still equals zero

    E = computeElementalMatrix(model, model.mets, false);
    S = model.S;
    [r,c]  = find(E'*S);
    if ~isempty(c)
        isValid = false;
        rxnsFound = model.rxns(unique(c));
        ignoreReactionsSel = findExcRxns(model);
        rxnSearchList = model.rxns(~ignoreReactionsSel);
        ignoreReactions = model.rxns(ignoreReactionsSel);
        for i=1:length(rxnSearchList)
            if strfind(rxnSearchList{i}, 'Ec_biomass')
                ignoreReactions = [ignoreReactions; rxnSearchList{i}];
            end
        end
        rxnsFound(ismember(rxnsFound, ignoreReactions)) = [];
        if isempty(rxnsFound)
            isValid = true;
        end
    else
        isValid = true;
        rxnsFound = [];
    end
end