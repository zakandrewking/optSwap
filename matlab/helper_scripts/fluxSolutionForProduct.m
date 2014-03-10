function [fluxes, soln, model] = fluxSolutionForProduct(model,target,swap)
    [model, biomass] = setupModel('iJO','EX_glc(e)','anaerobic','thko');
    if nargin > 2
        model = modelSwap(model, swap, false);
    end
    model.c(:) = 0;
    model.c(ismember(model.rxns,target)) = 1;
    model.lb(ismember(model.rxns,biomass)) = 0.1;
    soln = optimizeCbModel(model);
    fluxes = [model.rxns, num2cell(soln.x), num2cell(abs(soln.x))];
    fluxes = sortcell(fluxes, 3);
    fluxes = fluxes(:,1:2);
end
