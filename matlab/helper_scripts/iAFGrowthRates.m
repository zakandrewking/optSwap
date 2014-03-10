function iAFGrowthRates

% IAFGROWTHRATES
% 
% Calculate growth rates for iAF model, as
% shown in Figure 2 of Feist 2010.
% 
% 
% Zachary King 2012
    
    
    cd ~/Dropbox/lab/OptSwap/code
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% set cobra solver
    setenv 'GUROBI_HOME' '/Library/gurobi461/mac64'
    setenv 'PATH' '/Library/gurobi461/mac64/bin'
    setenv 'LD_LIBRARY_PATH' '/Library/gurobi461/mac64/lib:/usr/local/lib'
    setenv 'DYLD_LIBRARY_PATH' '/usr/local/lib'
    setenv 'GRB_LICENSE_FILE' '/Library/gurobi461/gurobi.lic'
    changeCobraSolver('glpk','LP');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%% load model
    model = loadModel('Ec_iAF1260_flux1');
    
   
    %%% glucose uptake to 20
    model = changeRxnBounds(model, 'EX_glc(e)', 0, 'l');
    
    %%% anaerobic
    model = changeRxnBounds(model, 'EX_o2(e)', 0, 'l');

    subs = {'EX_glc(e)';'EX_xyl_D(e)';'EX_glyc(e)'};
    subNames = model.rxnNames(ismember(model.rxns,subs));
    
    anGRs = cell(length(subs),1);
    for i = 1:length(subs) 
        modelTemp = changeRxnBounds(model, subs{i}, -20, 'l'); 
        soln = optimizeCbModel(modelTemp);
        anGRs{i} = soln.f;
    end
    display('anaerobic');
    display([subs anGRs]);
    
    
    
    %%% aerobic
    model = changeRxnBounds(model, 'EX_o2(e)', -20, 'l');
    
    aeGRs = cell(length(subs),1);
    for i = 1:length(subs) 
        modelTemp = changeRxnBounds(model, subs{i}, -20, 'l'); 
        soln = optimizeCbModel(modelTemp);
        aeGRs{i} = soln.f;
    end
    display('aerobic');
    display([subs aeGRs]);
    
    save_to_base(1);
    
end
