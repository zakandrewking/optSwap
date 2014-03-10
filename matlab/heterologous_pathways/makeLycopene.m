function model = makeLycopene(model)
% add reactions for lycopene production
% crtE, crtI, crtB genes from E. herbicola in pK19, KmR

    if false
       model = addReaction(model, 'lyc_rxn', ...
                 {'g3p[c]','pyr[c]','nadph[c]','ctp[c]','atp[c]',...
                  'co2[c]','nadp[c]','cmp[c]','adp[c]','ppi[c]'},...
                 [-8, -8, -16, -8, -8, 8, 16, 8, 8, 12], 0);
                 model.c(model.c~=0) = 0;
                 model.c(ismember(model.rxns,'lyc_rxn')) = 1;
                 return
    end
    
    newNames = {% 'FPS',
               'crtE',
               'crtB',
               'crtI',
                'lycotex',
               'EX_lyco(e)'};
    newMets = {% {'ipdp[c]','ppi[c]', 'grdp[c]'},
               {'ipdp[c]','frdp[c]','ggpp[c]','ppi[c]'},
               {'ggpp[c]', 'phyto[c]', 'ppi[c]'},
               {'phyto[c]', 'nadp[c]', 'lyco[c]', 'nadph[c]'},
              {'lyco[c]','lyco[e]'},
              {'lyco[e]'}};
    newRxns = {% [-2, 1, 1],
               [-1, -1, 1, 1],
               [-2, 1, 2],
               [-1, -8, 1, 8],
               [-1, 1],
               [-1]};
    rev = [% 0,
           0, 0, 0, 0, 0];
    
    for i=1:length(newNames) 
        model = addReaction(model, newNames{i}, newMets{i}, ...
                            newRxns{i}, rev(i));
        if i==2
            model.metFormulas(ismember(model.mets,'ggpp[c]')) = {'C20H33O7P2'};
        elseif i==3
            model.metFormulas(ismember(model.mets,'phyto[c]')) = {'C40H64'};
        elseif i==4
            model.metFormulas(ismember(model.mets,'lyco[c]')) = {'C40H56'};
        elseif i==5
            model.metFormulas(ismember(model.mets,'lyco[e]')) = {'C40H56'};
        end
            
        [isValid, foundRxns] = checkESMatrix(model);
        display(foundRxns)
    end
    
    % model.c(model.c~=0) = 0;
    % model.c(ismember(model.rxns,'EX_lyco(e)')) = 1;
end

% '4c2me[c]' '4-(cytidine 5'-diphospho)-2-C-methyl-D-erythritol' / 4-diphophocytidyl-2-C-methyl-d-erythritol
% '2p4c2me[c] 2-phospho-4-(cytidine 5'-diphospho)-2-C-methyl-D-erythritol / 4-diphosphocytidyl-2-C-methyl-2-phosphate-d-erythritol
% 'dmpp[c]' Dimethylallyl diphosphate/Dimethylallyl pyrophosphate
% 'dxyl5p[c]' '1-deoxy-D-xylulose 5-phosphate'
% 'frdp[c]' Farnesyl diphosphate/trans, trans Farnesyl pyrophosphate 
% 'grdp[c]' trans Geranyl pyrophosphate
% 'ggpp[c]' GGPP: Geranylgeranyl PP [new] C20H36O7P2
% 'h2mb4p[c]' '1-hydroxy-2-methyl-2-(E)-butenyl 4-diphosphate'
% 'ipdp[c]' 'Isopentenyl diphosphate'
% 'lyco[c]' LYCO: 'Lycopene' [new] C40H56
% '2mecdp[c]' '2-C-methyl-D-erythritol 2,4-cyclodiphosphate'
% '2me4p[c]' '2-C-methyl-D-erythritol 4-phosphate'
% 'octdp[c]' 'all-trans-Octaprenyl diphosphate'
% 'phyto[c]' PHYTO: Phytoene [new] C40H64
% 'ppi[c]' disphosphate/pyrophosphate
% 'pyr[c]' pyruvate
% 'g3p[c]' 'Glyceraldehyde 3-phosphate'
% 'udcpdp[c|p]' 'Undecaprenyl diphosphate'
% 
% is there evidence for this reaction?:
% Farnesyl pyrophosphate synthetase / ispA / 2IPPP → GPP + PPI
% (Alper, 2005)
% 
% rxns:
% DXPS '1-deoxy-D-xylulose 5-phosphate synthase'
% IPDDI 'isopentenyl-diphosphate D-isomerase'
% 'DMPPS' '1-hydroxy-2-methyl-2-(E)-butenyl 4-diphosphate reductase (dmpp)'
% 'h[c] + h2mb4p[c] + nadh[c]  -> dmpp[c] + h2o[c] + nad[c] '
% 'IPDPS' '1-hydroxy-2-methyl-2-(E)-butenyl 4-diphosphate reductase (ipdp)'
% 'h[c] + h2mb4p[c] + nadh[c]  -> h2o[c] + ipdp[c] + nad[c] '