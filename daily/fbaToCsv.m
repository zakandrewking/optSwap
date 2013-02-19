substrate = 'EX_glc(e)';
aerobicString = 'aerobic';
target = 'EX_ala-L(e)';
filename = '3-swap-ala-design.csv';
kos = {'ACKr';'ALCD2x';'LDH_D';'RPE'};
swaps = [];

model = setupModel('iJO',substrate,aerobicString,'nothko');
model = changeRxnBounds(model, kos, 0, 'b');
if ~isempty(swaps)
    model = modelSwap(model,swaps,false);
end
soln = optimizeCbModel(model,'max',0,0);
output = [model.rxns num2cell(soln.x)];
cell2csv(filename,output);
save_to_base(1)