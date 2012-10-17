function loopOptSwapDeded

% setup Cleaner
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = '1 swap on all max yield results about 10%';

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
                  'EX_ade(e)';};
    sets = [0,1;];
    for i=1:size(sets,1)
        status = sprintf('run: %d, kos: %d swaps: %d', i, sets(i,1), sets(i,2));
        opt.knockoutNum = sets(i,1);
        opt.swapNum = sets(i,2);
        opt.targetRxns = targetRxns;
        opt.experiment = run;
        opt.logFile = 'database-3.csv';
        opt.swapAllDhs = false;
        runOptSwapD(opt);
    end
    status = 'finished';
end