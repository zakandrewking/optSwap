function [newModel, newNames, coupling] = modelSwap(model,dhRxns,keepWtRxn)

% MODELSWAP
% 
% Return a model with swapped reactions
% 
% INPUTS
% model
% dhRxns
% keepWtRxn
% 
% OUTPUTS
% newModel
% newNames
% coupling
% 
% Zachary King 9/6/2012
    
    if isempty(dhRxns)
        newModel = model; newNames = []; coupling = [];
        return
    end
    if ~iscell(dhRxns)
        dhRxns = {dhRxns};
    end

    newNames = cell(size(dhRxns));
    
    for i=1:length(dhRxns)

        % turn off the old dhRxn
        if ~keepWtRxn
            model = changeRxnBounds(model, dhRxns{i}, 0, 'b');
            display('Setting old reaction bounds to zero');
        else
            display('Keeping old reaction bounds');
        end
        
        S = full(model.S);
        rxn = S(:,find(ismember(model.rxns,dhRxns{i})));
        nadInd = find(ismember(model.mets,'nad[c]'));
        nadhInd = find(ismember(model.mets,'nadh[c]'));
        nadpInd = find(ismember(model.mets,'nadp[c]'));
        nadphInd = find(ismember(model.mets,'nadph[c]'));
        newRxn = rxn;

        if (rxn(nadInd))
            newRxn(nadpInd) = rxn(nadInd);
            newRxn(nadInd) = 0;
        end
        if (rxn(nadhInd))
            newRxn(nadphInd) = rxn(nadhInd);
            newRxn(nadhInd) = 0;
        end
        if (rxn(nadpInd))
            newRxn(nadInd) = rxn(nadpInd);
            newRxn(nadpInd) = 0;
        end
        if (rxn(nadphInd))
            newRxn(nadhInd) = rxn(nadphInd);
            newRxn(nadphInd) = 0;
        end

        newMets = model.mets(newRxn~=0);
        newRxn = newRxn(newRxn~=0);
        rev = model.rev(ismember(model.rxns,dhRxns{i}));
        newNames{i} = [dhRxns{i} 'swap'];

        % check if reaction already exists
        selRxn = zeros(length(model.mets),1);
        for j=1:length(newMets)
            selRxn(find(ismember(model.mets,newMets{j})),1) = newRxn(j);
        end
        indOfExisting = find(ismember(model.S',selRxn','rows'));
        if ~isempty(indOfExisting)
            if length(indOfExisting)>1, error('multiple copies of reaction found!'); end
            existingName = model.rxns(indOfExisting(1));
            disp(sprintf('Swap of reaction %s already exists: %s',dhRxns{i}, ...
                 existingName{1}));
            newNames{i} = existingName{1};
        else
            model = addReaction(model, newNames{i},...
                            newMets,newRxn,rev);
        end
        
        coupling(i,1:2) = [find(ismember(model.rxns, dhRxns{i})),...
                                               find(ismember(model.rxns, newNames{i}))];
    end
    
    % check for double ups
    a = intersect(coupling(:,1), coupling(:,2));
    if ~isempty(a)
        error(sprintf(['found %d double entries. For clarity, remove swaps from ' ...
                       'list of dehydrogenase reactions.'], length(a)/2));
    end
    newModel = model;

    [isValid, foundRxns] = checkESMatrix(newModel);
    if (~isValid)
        display(foundRxns)
        warning(sprintf('invalid S matrix. E*S ~= 0. Found %d bad rxns.', ...
                      length(foundRxns)))
    end
end
