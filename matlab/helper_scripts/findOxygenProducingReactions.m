function [reactions, formula, flux] = findOxygenProducingReactions(model)
% give a model, find all reactions that are producing oxygen

    % run pFBA
    [~,~,modelIrrevFM] = pFBA(modelTemp,'geneoption',0,'tol',1e-7);
    soln = optimizeCbModel(modelIrrevFM);

    % get reversibility
    [modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(model);

    % find o2 producing reactions
    o2_met = 's_1275'; % oxygen [cytoplasm]
    modelIrrevFM.S(ismember(model.mets, o2_met), :);
end