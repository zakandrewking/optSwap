function singleCarrier

model = loadModelNamed('iJO');
% anaerobic
model = changeRxnBounds(model, 'EX_o2(e)', 0, 'l');

% % minimum biomass production
% model.lb(model.c~=0) = 0.1;

modelWT = model;

% remove transhydrogenases
thRxns = {'NADTRHD', 'THD2pp'};
model = changeRxnBounds(model, thRxns, 0, 'b');


model = replaceNadpWithNad(model);

soln   = optimizeCbModel(model);
solnWT = optimizeCbModel(modelWT);

biomassRxn = model.rxns(model.c~=0);
figure(1)
productionEnvelope(modelWT, [], '-k', 'EX_for(e)', biomassRxn);
figure(2)
productionEnvelope(model, [], '-k', 'EX_for(e)', biomassRxn);

display(sprintf('wild type biomass:        %f', solnWT.f));
display(sprintf('single e carrier biomass: %f', soln.f));


end


function model = replaceNadpWithNad(model)

    nadph_ind = find(ismember(model.mets,'nadph[c]'));
    nadp_ind  = find(ismember(model.mets,'nadp[c]'));
    nadh_ind  = find(ismember(model.mets,'nadh[c]'));
    nad_ind   = find(ismember(model.mets,'nad[c]'));
    
    model.S(nadh_ind,:) = model.S(nadh_ind,:) + model.S(nadph_ind,:);
    model.S(nad_ind,:) = model.S(nad_ind,:) + model.S(nadp_ind,:);
    model = removeMetabolites(model, model.mets([nadph_ind, nadp_ind]));

end