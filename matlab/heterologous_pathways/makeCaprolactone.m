function model = makeCaprolactone(model)
% add reactions for lycopene production
% crtE, crtI, crtB genes from E. herbicola in pK19, KmR

    if true
       model = addReaction(model, 'CMHO', ...
                 {'nadph[c]','h[c]','nadp[c]'},...
                 [-1 1 1], 0);
                 model.c(model.c~=0) = 0;
                 model.c(ismember(model.rxns,'CMHO')) = 1;
    end

end