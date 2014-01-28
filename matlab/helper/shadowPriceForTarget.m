function [mets, prices] = shadowPriceForTarget(model, biomass, target, swap)

    mets = {'glc-D[c]'; 'nadp[c]'; 'nadph[c]'; 'h[c]'; 'for[c]'};
    
    if nargin > 3
        model = modelSwap(model, swap, false);
    end
    model.c(:) = 0;
    model.c(ismember(model.rxns,target)) = 1;
    model.lb(ismember(model.rxns,biomass)) = 0.1;
    soln = optimizeCbModel(model);
    
    sel = ismember(model.mets, mets);
    mets = model.mets(sel);
    prices = soln.w(sel);
    
end