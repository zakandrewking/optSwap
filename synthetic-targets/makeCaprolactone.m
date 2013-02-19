function model = makeCaprolactone(model)
% add reactions for lycopene production
% crtE, crtI, crtB genes from E. herbicola in pK19, KmR

    if true
       model = addReaction(model, 'cap_rxn', ...
                 {'o2[c]','nadph[c]','h[c]',...
                  'h2o[c]','nadp[c]'},...
                 [-1 -2 -2 2 2], 0);
                 model.c(model.c~=0) = 0;
                 model.c(ismember(model.rxns,'cap_rxn')) = 1;
    end

end