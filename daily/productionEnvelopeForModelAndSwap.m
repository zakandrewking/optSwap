figure()
aerobicString = 'aerobic';
[model_thko,biomass] = setupModel('iJO','EX_glc(e)',aerobicString,'thko');
model_wt = setupModel('iJO','EX_glc(e)',aerobicString,'nothko');
model_wt = makeCaprolactone(model_wt);
model_thko = makeCaprolactone(model_thko);
model_wt_gapd = modelSwap(model_wt,'GAPD',false);
model_thko_gapd = modelSwap(model_thko,'GAPD',false);

productionEnvelope(model_thko_gapd,[],'b','cap_rxn',biomass);
hold on;
productionEnvelope(model_thko,[],'r','cap_rxn',biomass);
productionEnvelope(model_wt_gapd,[],'--b','cap_rxn',biomass);
productionEnvelope(model_wt,[],'--r','cap_rxn',biomass);
legend('gapd swap, thko','no swap, thko','gapd swap, no thko','no swap, no thko');
set(gcf,'Color','White');
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Units', 'pixels');
set(gcf, 'Position', [500 500 400 280]);
xlabel('growth rate (h^{-1})');
ylabel('\epsilon-caprolactone production (mmol gDW^{-1} h^{-1})');
title([aerobicString ' production of \epsilon-caprolactone']);
