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
    targetRxns = {'EX_etoh(e)';
                  'EX_lac-D(e)';
                  'EX_ala-L(e)';};
    % turn on special cases
    % model = makeFumTransporterReversible(model);

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


