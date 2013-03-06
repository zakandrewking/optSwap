function model = loadModelNamed(modelName)

% LOADMODEL
%
% loads SBML model and saves as .mat file for future use
%
% Zachary King 2012

    addRxns = false;
    c = pathsep;
    if (strcmp(modelName,'iAF') || strcmp(modelName, 'iAF1260'))
        modelName = 'Ec_iAF1260_flux1';
    elseif strcmp(modelName,'iJO')
        modelName = 'iJO1366';
    elseif strcmp(modelName,'iAF1260b')
        modelName = 'Ec_iAF1260_flux1';
        addRxns = true;
    elseif strcmp(modelName,'iND750')
        modelName = 'Sc_iND750_flux1';
    elseif strcmp(modelName,'iJO-h')
        modelName = 'iJO1366-heterogenous-pathways';
    end


    try
        model = load([modelName '.mat']);
    catch
        model = readCbModel(['/Users/zaking/Dropbox/git/Matlab/models/' modelName],1000,'SBML'); 
        save([modelName '.mat'],'-struct','model');
    end

    if addRxns
        % remove the old rxn
        model = removeRxns(model,{'GLYt2pp', 'ASPt2pp', 'ALAt2pp'});

        % add reactions from Feist 2010 Supp Table 1...iAF1260b
        model = addReaction(model, ...
                            'DHORDfum', 'dhor_S[c] + fum[c]  -> orot[c] + succ[c] ');
        % is this DHORTfum from Feist 2010?
        model = addReaction(model, ...
                            'MALt3pp', 'mal_L[c] + h[p]  -> h[c] + mal_L[p]');
        model = addReaction(model, ...
                            'ALAt2rpp', 'ala_L[p] + h[p]  <=> ala_L[c] + h[c]');
        model = addReaction(model, ...
                            'GLYt2rpp', 'gly[p] + h[p]  <=> gly[c] + h[c]');
        model = addReaction(model, ...
                            'CITt3pp', 'cit[c] + h[p]  -> h[c] + cit[p]');
        model = addReaction(model, ...
                            'ASPt2rpp', 'asp_L[p] + h[p]  <=> asp_L[c] + h[c]');
    end
end