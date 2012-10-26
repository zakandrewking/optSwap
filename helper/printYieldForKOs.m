function printYieldForKOs(targets,kos,substrate,aerobicString)

    model = setupModel('iJO',substrate,aerobicString,'thko');
    out = zeros(length(targets));
    for i=1:length(targets)
        rxn = targets{i};
        modelT = model;
        modelT = setupModelForTarget(modelT, rxn);
        modelT = changeRxnBounds(modelT, kos{i}, 0, 'b');
        soln = optimizeCbModel(modelT);
        if ~isempty(soln.x)
            out(i) = soln.x(ismember(modelT.rxns, rxn));
        end
    end 
    for i =1:length(targets) 
        fprintf('%.4f\n', out(i));
    end
end