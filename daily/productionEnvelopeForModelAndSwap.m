 subSystems = {' 1,3-Propanediol production';
                  ' 1,4-Butanediol production';
                  ' 2,3-Butanediol production';
                  ' 2,3-Butanediol production';
                  ' 2,3-Butanediol production';
                  ' 3-Hydroxybutyrate production';
                  ' 3-Hydroxybutyrate production';
                  ' 3-Hydroxypropanoate production';
                  ' 3-Hydroxyvalerate production';
                  ' 3-Hydroxyvalerate production';
                  'Styrene production';
                  'Styrene production';
                  'p-Hydroxystyrene production'};

figure()
hold on;
aerobicString = 'anaerobic';
[model_thko,biomass] = setupModel('iJO-h','EX_glc(e)',aerobicString,'thko');
model_wt = setupModel('iJO-h','EX_glc(e)',aerobicString,'nothko');
model_thko = turnOnSubSystem(model_thko, 8, subSystems);
model_wt = turnOnSubSystem(model_wt, 8, subSystems);
model_thko_gapd = modelSwap(model_thko,'GAPD',false);
model_thko_gapd_pdh = modelSwap(model_thko,{'GAPD','PDH'},false);

productionEnvelope(model_thko_gapd,[],'b','EX_3hpp(e)',biomass);
productionEnvelope(model_thko_gapd_pdh,[],'r','EX_3hpp(e)',biomass);
productionEnvelope(model_thko,[],'k','EX_3hpp(e)',biomass);
productionEnvelope(model_wt,[],'--k','EX_3hpp(e)',biomass);
legend('gapd swap, thko','gapd pdh swap, thko','thko','wt');
set(gcf,'Color','White');
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Units', 'pixels');
set(gcf, 'Position', [500 500 400 280]);
xlabel('growth rate (h^{-1})');
