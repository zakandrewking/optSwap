function fermenationShift
    dhRxns = {{'GAPD','G6PDH2r'};};
    model = setupModel('iJO','EX_glc(e)','anaerobic','noTHKO');
    out = zeros(length(dhRxns),6);
    for i=1:length(dhRxns)
        modelT = model;
        if ~isempty(dhRxns{i})
            %modelT = modelSwap(modelT, dhRxns{i}, false);
            modelT = changeRxnBounds(modelT, dhRxns{i}, 0, 'b');
        end
        soln = optimizeCbModel(modelT);
        
        sel = ismember(modelT.rxns,{'EX_for(e)', 'EX_ac(e)', 'EX_etoh(e)',...
                            'EX_succ(e)', 'EX_akg(e)', 'EX_lac-D(e)'});
        if ~isempty(soln.x)
            a = soln.x(sel)';
        else
            a = zeros(1,6);
        end
        out(i,:) = a;          
    end
    disp(out)
end