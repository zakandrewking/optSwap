function model = makeIsobutanol(model)
% make isobutanol 
% 1. Trinh CT. Elucidating and reprogramming Escherichia coli
% metabolisms for obligate anaerobic n-butanol and isobutanol
% production. Applied microbiology and biotechnology. 2012;95(4):1083–94.
% 
% IBUT1 (acetolactate synthase, AlsSBS): 2 pyruvate -> acetolactate + CO2 
%    = ACLS in ijo1366
% IBUT2 (2,3-dihydroxy isovalerate oxidore- ductase, IlvCEC): acetolactate+NADPH+H+
% 02,3- dihydroxyvalerate+NADP+; 
% IBUT3 (2,3-dihydroxy isoval- erate dehydratase,
% IlvDEC): 2,3-dihydroxy isovalerate0iso- valerate+H2O; 
% IBUT4 (ketoacid
% decarboxylase, KivdLL): isovalerate0isobutanal+CO2; and 
% IBUT5 (isobutanal de-
% hydrogenase, AdhE2CA): isobutanal+NADH+H+0isobu- tanol+NAD+

    newNames = {'IBUT2',
                'IBUT3',
                'IBUT4',
                'IBUT5',
                'EX_ibut(e)'};
    newMets = {{'alac-S[c]','h[c]','nadph[c]','dhv[c]','nadp[c]'},
               {'dhv[c]','ival[c]','h2o[c]'},
               {'ival[c]','h[c]','ibuta[c]','co2[c]'},
               {'ibuta[c]','nadh[c]','h[c]','ibut[c]','nad[c]'},
               {'ibut[c]'}};
    newRxns = {[-1, -1, -1, 1, 1],
               [-1, 1, 1],
               [-1, -1, 1, 1],
               [-1, -1, -1, 1, 1],
               [-1]};
    rev = [0, 0, 0, 0, 0];
    
    for i=1:length(newNames) 
        model = addReaction(model, newNames{i}, newMets{i}, ...
                            newRxns{i}, rev(i));
        if i==1
            model.metFormulas(ismember(model.mets,'dhv[c]')) = {'C5H9O4'};
        elseif i==2
            model.metFormulas(ismember(model.mets,'ival[c]')) = {'C5H7O3'};
        elseif i==3
            model.metFormulas(ismember(model.mets,'ibuta[c]')) = {'C4H8O'};
        elseif i==4
            model.metFormulas(ismember(model.mets,'ibut[c]')) = {'C4H10O'};
        end 
    end
    [isValid, foundRxns] = checkESMatrix(model)
end

% acetolactate ... 'alac-S[c]' ... '(S)-2-Acetolactate' ... C5H7O4
% 2,3-dihydroxyvalerate ... dhv[c] ... 2,3-dihydroxy-3-methylbutanoate ... C5H9O4
% isovalerate ...ival[c] ... 3-methyl-2-oxobutanoate ...  C5H7O3
% isobutanal ... ibuta[c] ... C4H8O
% isobutanol ... ibut[c] ... C4H10O
