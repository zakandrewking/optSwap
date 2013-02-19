function model = makePhb(model)
% make poly(3-hydroxybutyric acid) (PHB)
% 
% phaAB from Wautersia eutropha (formerly Ralstonia eutropha)
% phaEC from Allochromatium vinosm (formerly Chromatium Vinosum)
% (i) 2 AcCoA→AcAcCoA+CoASH; 
% (ii) AcAcCoA+NADPH→3HBCoA+NADP; and 
% (iii) 3HBCoA→PHB+CoASH.
% 
% 1. Tyo KEJ, Fischer CR, Simeon F, Stephanopoulos G. Analysis of
% polyhydroxybutyrate flux limitations by systematic genetic and metabolic
% perturbations. Metabolic engineering. 2010;12(3):187–95.
    
    
    
    newNames = {'phaB', 'phaC', 'EX_phb(e)'};
    newMets = {% {'accoa[c]','aacoa[c]','coa[c]'},
               {'aacoa[c]','nadph[c]','h[c]','3hbcoa[c]','nadp[c]'},
               {'3hbcoa[c]','phb[c]','coa[c]'},
               {'phb[c]'}};
    newRxns = {% [-2, 1, 1],
               [-1, -1, -1, 1, 1],
               [-1, 1, 1],
               [-1]};
    rev = [0,0,0]

    for i=1:length(newNames)
        model = addReaction(model, newNames{i}, newMets{i}, ...
                            newRxns{i}, rev(i));
        if i==2
            model.metFormulas(ismember(model.mets,'3hbcoa[c]')) = {'C25H38N7O18P3S'};
        elseif i==3
            model.metFormulas(ismember(model.mets,'phb[c]')) = {'C4H6O2'};
        end

        [isValid, foundRxns] = checkESMatrix(model);
        display(foundRxns)
    end

end

% 'gamma-hydroxybutyrate' == 4-hydroxybutanoate
% 3hbcoa[c] ... (R)-3-hydroxybutanoyl-CoA
% phb[c] ... poly-β-hydroxybutyrate
% 
% Warning: Model already has the same reaction you tried to add: ACACT1r
% 
% 'aacoa[c]'     'C25 H36 N7 O18 P3 S'
% '3hbcoa[c]'    'C25 H38 N7 O18 P3 S'
% 'phb[c]'       'C4   H6     O2'        
% 'coa[c]'       'C21 H32 N7 O16 P3 S'