function maxYield
% maxYield
% calculate maximum yield of target reactions at minimum biomass production
% with dehydrogenase swaps
%
% Zachary King 9/20/2012

    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'max-yield_all-forward-transporters_w-NADTRHD-open-case';

    logFile = sprintf('%s-%s.csv', run, datestr(now, 'yy-mm-dd_HH_MM_SS'));
    fileId = fopen(logFile, 'w');
    fprintf(fileId, '%s\n', logFile);
    fprintf(fileId, 'target, isAerobic, substrate, ko swap case\n')

    % set parameters
    dhRxns = dhRxnList(31);
    growthMin = 0.1;
    substrateList = {'EX_glc(e)', 'EX_xyl-D(e)'};
    isAerobic = [0,1];
    thko_swap = {'wt','thko','NADTRHD open'};

    % load model
    [model, biomassRxn] = setupModel('iJO', 'none', 'aerobic', 'noTHKO');
    model = changeRxnBounds(model, biomassRxn, growthMin, 'l');

    % get all exchange reactions
    targetRxns = {'EX_h2s(e)';
                  'EX_cys-L(e)';
                  'EX_hom-L(e)';
                  'EX_thr-L(e)';
                  'EX_15dap(e)';
                  'EX_lys-L(e)';
                  'EX_ile-L(e)';
                  'EX_leu-L(e)';
                  'EX_gly(e)';
                  'EX_5mtr(e)';
                  'EX_spmd(e)';
                  'EX_pro-L(e)';
                  'EX_cgly(e)';
                  'EX_asp-L(e)';
                  'EX_gthrd(e)';
                  'EX_ser-L(e)';
                  'EX_thym(e)';
                  'EX_ptrc(e)';
                  'EX_glyc(e)';
                  'EX_pyr(e)';
                  'EX_orn(e)';
                  'EX_agm(e)';
                  'EX_urea(e)';
                  'EX_arg-L(e)';
                  'EX_glyald(e)';
                  'EX_lac-L(e)';
                  'EX_h2o(e)';
                  'EX_asn-L(e)';
                  'EX_thymd(e)';
                  'EX_LalaDgluMdapDala(e)';
                  'EX_acser(e)';
                  'EX_ac(e)';
                  'EX_LalaDgluMdap(e)';
                  'EX_phe-L(e)';
                  'EX_tyr-L(e)';
                  'EX_ala-D(e)';
                  'EX_alaala(e)';
                  'EX_12ppd-S(e)';
                  'EX_his-L(e)';
                  'EX_g3pe(e)';
                  'EX_g3pg(e)';
                  'EX_enter(e)';
                  'EX_feenter(e)';
                  'EX_4abut(e)';
                  'EX_ura(e)';
                  'EX_ade(e)';
                  'DM_5DRIB';
                  'EX_pheme(e)';
                  'EX_adn(e)';
                  'EX_glyc3p(e)';
                  'EX_etha(e)';
                  'EX_glu-L(e)';
                  'EX_uri(e)';
                  'EX_indole(e)';
                  'EX_trp-L(e)';
                  'EX_cytd(e)';
                  'EX_gua(e)';
                  'EX_hxan(e)';
                  'DM_OXAM';
                  'EX_xan(e)';
                  'EX_alltn(e)';
                  'EX_acolipa(e)';
                  'EX_xtsn(e)';
                  'EX_colipa(e)';
                  'EX_ins(e)';
                  'EX_kdo2lipid4(e)';
                  'EX_acald(e)';
                  'EX_enlipa(e)';
                  'EX_anhgm(e)';
                  'EX_eca4colipa(e)';
                  'EX_lipa_cold(e)';
                  'EX_lipa(e)';
                  'EX_colipap(e)';
                  'EX_glyc-R(e)';
                  'EX_cit(e)';
                  'EX_akg(e)';
                  'EX_h2(e)';
                  'EX_glyclt(e)';
                  'EX_mal-L(e)';
                  'EX_fe2(e)';
                  'DM_MTHTHF';
                  'EX_fum(e)';
                  'EX_quin(e)';
                  'EX_dha(e)';
                  'EX_h(e)';
                  'EX_co2(e)';
                  'EX_12ppd-R(e)';
                  'EX_succ(e)';
                  'EX_for(e)';};
    % turn on special cases
    model = makeFumTransporterReversible(model);

    for m=1:length(targetRxns)
        model_m = model;
        model_m.c = zeros(size(model_m.c));
        model_m.c(ismember(model_m.rxns,targetRxns{m})) = 1;
        for i=1:length(isAerobic)
            model_a = model_m;
            if isAerobic(i)==0, model_a=changeRxnBounds(model_a,'EX_o2(e)',0,'b'); end
            for j=1:length(substrateList)
                model_s = model_a;
                model_s = changeRxnBounds(model_s, substrateList{j}, -20, 'l');
                for k=1:length(thko_swap)
                    model_t = model_s;
                    switch thko_swap{k}
                      case 'thko'
                        model_t=knockoutTH(model_t);
                      case 'NADTRHD open'
                        model_t=changeRxnBounds(model_t,'NADTRHD',-1000,'l');
                        model_t=changeRxnBounds(model_t,'NADTRHD',1000,'u');
                    end
                    soln = optimizeCbModel(model_t);
                    fprintf(fileId,'%s,%d,%s,%s,%.4f\n',targetRxns{m},isAerobic(i), ...
                            substrateList{j}, thko_swap{k}, soln.f);
                end
            end
        end
    end
    fclose(fileId);
end

function model = makeFumTransporterReversible(model)
    transporter = {'FUMt2_2pp'};
    model.rev(ismember(model.rxns,transporter)) = 1;
    model.lb(ismember(model.rxns,transporter)) = -1000;
end

function model = knockoutTH(model)
    thRxns = {'NADTRHD', 'THD2pp'};
    model = changeRxnBounds(model, thRxns, 0, 'b');
end


