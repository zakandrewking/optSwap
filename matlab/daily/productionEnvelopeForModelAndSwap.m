 % subSystems = {' 1,3-Propanediol production';
 %                  ' 1,4-Butanediol production';
 %                  ' 2,3-Butanediol production';
 %                  ' 2,3-Butanediol production';
 %                  ' 2,3-Butanediol production';
 %                  ' 3-Hydroxybutyrate production';
 %                  ' 3-Hydroxybutyrate production';
 %                  ' 3-Hydroxypropanoate production';
 %                  ' 3-Hydroxyvalerate production';
 %                  ' 3-Hydroxyvalerate production';
 %                  'Styrene production';
 %                  'Styrene production';
 %                  'p-Hydroxystyrene production'};

figure()
hold on;
aerobicString = 'anaerobic';
% [model_thko,biomass] = setupModel('iND750','EX_xyl-D(e)',aerobicString, ...
%                                            'nothko');
[model_wt, biomass] = setupModel('iND750','EX_xyl_D(e)',aerobicString,'nothko');
% model_thko = turnOnSubSystem(model_thko, 8, subSystems);
% model_wt = turnOnSubSystem(model_wt, 8, subSystems);
% model_thko_gapd = modelSwap(model_thko,'GAPD',false);
model_no_thko_gapd = modelSwap(model_wt, 'GAPD', false);
% model_thko_gapd_pdh = modelSwap(model_thko,{'GAPD','PDH'},false);

% productionEnvelope(model_thko_gapd,[],'b','EX_cys-L(e)',biomass);
% productionEnvelope(model_thko_gapd_pdh,[],'r','EX_cys-L(e)',biomass);
productionEnvelope(model_no_thko_gapd, [], 'r', 'EX_cys_L(e)', biomass);
% productionEnvelope(model_thko,[],'k','EX_cys-L(e)',biomass);
productionEnvelope(model_wt,[],'k','EX_cys_L(e)',biomass);
legend('gapd swap, nothko','wt');
set(gcf,'Color','White');
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Units', 'pixels');
set(gcf, 'Position', [500 500 400 280]);
xlabel('growth rate (h^{-1})');
ylabel('L-Cysteine production');
filename = 'iND750--Cys production on D-xylose';
title(filename);
print('-depsc2', filename);