figure()
clear
[model,biomass] = setupModel('iJO','EX_glc(e)','anaerobic');

targets = {'lyc_rxn';
           'cap_rxn';
           'EX_phb(e)';
           'EX_ibut(e)'};

% targets = {'EX_cys-L(e)';
%            'EX_hom-L(e)';
%            'EX_thr-L(e)';
%            'EX_15dap(e)';
%            'EX_lys-L(e)';
%            'EX_val-L(e)'};
           
%            'EX_ile-L(e)';
%            'EX_leu-L(e)';
%            'EX_gly(e)';
%            'EX_5mtr(e)';
%            'EX_spmd(e)';
%            'EX_pro-L(e)';
%            'EX_cgly(e)';
%            'EX_asp-L(e)';
%            'EX_gthrd(e)';
%            'EX_ser-L(e)';
%            'EX_thym(e)';
%            'EX_ptrc(e)';
%            'EX_glyc(e)';
%            'EX_pyr(e)';
%            'EX_orn(e)';
%            'EX_agm(e)';
%            'EX_urea(e)';
%            'EX_arg-L(e)';
%            'EX_glyald(e)';
%            'EX_lac-L(e)';
%           };
set(gca, 'ColorOrder', hsv(length(targets)))
hold all
% markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
for i=1:length(targets)
    if i==1, modelT = makeLycopene(model);
    elseif i==2, modelT = makeCaprolactone(model);
    elseif i==3, modelT = makePhb(model);
    elseif i==4, modelT = makeIsobutanol(model);
    else, modelT = model; end
    % options.marker = markers{mod(i,numel(markers))+1}
    oxygenEnvelope(modelT,targets{i});  
end
legend(targets, 'Interpreter', 'None');