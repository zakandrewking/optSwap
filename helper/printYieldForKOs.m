function printYieldForKOs(targets,swaps,kos,substrate,aerobicString)
% make swaps empty [] for just knocks
    noSwaps = isempty(swaps);
    if noSwaps, swaps = cell(size(kos)); end
    
    model = setupModel('iJO',substrate,aerobicString,'thko');
    out = zeros(length(targets),4);
    for i=1:length(targets)
        rxn = targets{i};
        if ~noSwaps
            modelT = model;
            modelT = setupModelForTarget(modelT, rxn);
            modelT = changeRxnBounds(modelT, kos{i}, 0, 'b');
            modelT = modelSwap(modelT, swaps{i}, false);
            soln = optimizeCbModel(modelT);
            if ~isempty(soln.x)
                out(i,1) = soln.x(ismember(modelT.rxns, rxn));
                out(i,3) = out(i,1) * soln.f;
            end
        end
        modelT = model;
        modelT = setupModelForTarget(modelT, rxn);
        modelT = changeRxnBounds(modelT, [kos{i}; swaps{i}], 0, 'b');
        soln = optimizeCbModel(modelT);
        if ~isempty(soln.x)
            out(i,2) = soln.x(ismember(modelT.rxns, rxn));
            out(i,4) = out(i,2) * soln.f;
        end
    end 
    if noSwaps
        fprintf('knocks / SSP knocks\n');
        for i =1:length(targets) 
            fprintf('%.4f\t%.4f\n', out(i,2), out(i,4));
        end
    else
        fprintf('swaps and knocks / SSP swaps and knocks/ knocks / SSP knocks\n');
        for i =1:length(targets) 
            fprintf('%.4f\t%.4f\t%.4f\t%.4f\n', out(i,1), out(i,3), out(i,2), out(i,4));
        end
    end
end