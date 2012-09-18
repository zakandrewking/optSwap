function [model, biomassRxn] = setupModel(modelName,substrate,isAerobic,transhydrogenaseKnockout)
% setupModel
% 
% INPUTS
%   modelName
%   substrate
%   isAerobic
%   transhydrogenaseKnockout
% 
% OUTPUTS
%   model
%   biomassRxn
%
% Zachary King 9/12/12
    
    % load the model
    model = loadModelNamed(modelName);
    biomassRxn = model.rxns(model.c~=0);
    
    % prepare the model
    if transhydrogenaseKnockout
        % knock out th's
        thRxns = {'NADTRHD', 'THD2pp'};
        model = changeRxnBounds(model, thRxns, 0, 'b');
    end
    model = changeRxnBounds(model, 'EX_glc(e)', 0, 'l');
    model = changeRxnBounds(model, substrate, -20, 'l');
    if isAerobic
        model = changeRxnBounds(model, 'EX_o2(e)', -20, 'l');
    else
        model = changeRxnBounds(model, 'EX_o2(e)', 0, 'l'); 
        % BAD REACTIONS, DOWN BOY
        model = changeRxnBounds(model, {'CAT';'SPODM';'SPODMpp'}, [0;0;0], 'b');
        % can make oxygen...turn off for anaerobic simulations 
    end
    % Swap reactions, keeping wild type dehydrogenases
    % [~, newNames] = modelSwap(model, dhRxns, true);

end