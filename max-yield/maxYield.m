function maxYield

    run = sprintf('max-yield_rem-loops_glc-4-targets_%s', datestr(now, 'yy-mm-dd_HH_MM_SS'));
    logFile = [run '.csv'];
    global fileId dhRxns newNames targetRxns growthMins fluxes swapFluxes
    fileId = fopen(logFile, 'a');

    model = loadModelNamed('iJO');

    targetRxns = {% 'EX_etoh(e)'
                  % 'EX_lac-D(e)'
                  % 'EX_glyc(e)'
                  'EX_ala-L(e)'
                  'EX_pyr(e)'
                  % 'EX_fum(e)'
                  % 'EX_succ(e)'
                  % 'EX_akg(e)',
                  'EX_for(e)',
                  'EX_ac(e)',
                 };

    dhRxns = {'3OAR100'
              '3OAR40'
              '3OAR60'
              '3OAR80'
              'AGPR'
              'AKGDH'
              'ASAD'
              'DHDPRy'
              'EAR40x'
              'EAR60x'
              'EAR80x'
              'G6PDH2r'
              'GAPD'
              'GLUDy'
              'GND'
              'HPYRRx'
              'HSDy'
              'ICDHyr'
              'IPMD'
              'KARA1'
              'KARA2'
              'MDH'
              'MTHFD'
              'NADH16pp'
              'PDH'
              'PGCD'
              'SHK3Dr'
              'TRSARr'
              'ME1'
              'ME2'
              'LDH_D'
              'LCARR'};
    % model = makeIrreversible(model,dhRxns);

    biomassRxn = model.rxns(model.c~=0);
    substrateList = {'EX_glc(e)',
                     % 'EX_glyc(e)',
                     % 'EX_xyl-D(e)',
                    };

    growthMins = [0.0, 0.1, 0.3];

    model = makeFumTransporterReversible(model);

    isAerobic = false;
    if isAerobic
        model = changeRxnBounds(model, 'EX_o2(e)', -20, 'l');
    else
        model = changeRxnBounds(model, 'EX_o2(e)', 0, 'l');
    end

    out = zeros(length(substrateList), length(targetRxns),...
                length(growthMins), 3);
    fluxes = cell(size(out)); swapFluxes = fluxes;

    for m = 3:3

        model_K = model;

        switch m
          case 1
            disp('1');
          case 2
            disp('2');
            % knock out th's
            thRxns = {'NADTRHD', 'THD2pp'};
            model_K = changeRxnBounds(model_K, thRxns, 0, 'b');
          case 3
            disp('3');
            thRxns = {'NADTRHD', 'THD2pp'};
            model_K = changeRxnBounds(model_K, thRxns, 0, 'b');
            [model_K, newNames] = modelSwap(model_K, dhRxns, true);
            % model_K = makeIrreversible(model_K, newNames);
        end

        for i=1:length(substrateList)

            substrate = substrateList{i};

            model_KS = model_K;
            model_KS = changeRxnBounds(model_KS, 'EX_glc(e)', 0, 'l');
            model_KS = changeRxnBounds(model_KS, substrate, -20, 'l');

            for j = 1:length(targetRxns)

                model_KST = model_KS;
                model_KST.c = zeros(size(model_KST.c));
                model_KST.c(ismember(model_KST.rxns, targetRxns{j})) = 1;

                for k = 1:length(growthMins)

                    model_KSTG = model_KST;

                    model_KSTG.lb(ismember(model_KSTG.rxns, biomassRxn)) = growthMins(k);


                    disp('running optimization with no loops');
                    changeCobraSolverParams('MILP', 'printLevel', 3); 
                    sol = optimizeCbModel(model_KSTG,'max',0,0);
                    if ~isempty(sol.x)
                        fluxes{i,j,k,m} = sol.x(ismember(model_KSTG.rxns,dhRxns));
                        if m==3
                            swapFluxes{i,j,k,m} = sol.x(ismember(model_KSTG.rxns, ...
                                                                 newNames));
                        end
                    else
                        fluxes{i,j,k,m} = [];
                    end

                    out(i,j,k,m) = sol.f;


                end
            end
        end

    end
    
    for i=1:length(substrateList)
        fprintf(fileId, '%s\n\n', substrate)
        fprintf(fileId, ',0.0,,,0.1,,,0.2,,,0.3\n,');
        for x=1:4
            fprintf(fileId, 'noThKO,THKO,THKO+swaps,');
        end
        fprintf(fileId, '\n');
        for j=1:length(targetRxns)
            fprintf(fileId, '%s,', targetRxns{j});
            for k=1:length(growthMins)
                for m=1:3
                    fprintf(fileId, '%.2f,', out(i,j,k,m));
                end
            end
        end
    end
    
    fprintf(fileId, '\n\n\nFLUXES\n\n');
      
    for i=1:length(substrateList)
        printSubstrate(substrateList{i},i);
    end

    fclose(fileId);

end


function printSubstrate(substrate,i)
    global fileId targetRxns
    fprintf(fileId, '%s\n\n', substrate);

    for j=1:length(targetRxns)
        printTarget(targetRxns{j},i,j);
    end

    fprintf(fileId, '\n');
end

function printTarget(target,i,j)
    global fileId dhRxns newNames
    fprintf(fileId, 'target: %s\n', target);
    fprintf(fileId, ',0.0,,,0.1,,,0.2,,,0.3\n,');
    for x=1:4
        fprintf(fileId, 'noThKO,THKO,THKO+swaps,');
    end
    fprintf(fileId, '\n');
    for n=1:length(dhRxns)
        printDhRxn(dhRxns{n},i,j,n,false);
    end
    for n=1:length(newNames)
        printDhRxn(newNames{n},i,j,n,true);
    end
end

function printDhRxn(dhRxn,i,j,n,swap)
    global fileId growthMins fluxes swapFluxes
    fprintf(fileId, '%s,', dhRxn);
    for k=1:length(growthMins)
        for m=1:3
            if swap && (m < 3)
                fprintf(fileId, ',');
                continue;
            elseif swap
                fl = swapFluxes{i,j,k,m};
            else
                fl = fluxes{i,j,k,m};
            end
            if isempty(fl), fl = 0;
            else fl = fl(n); end
            fprintf(fileId, '%.2f,', fl);
        end
    end
    fprintf(fileId, '\n');
end


function model = makeFumTransporterReversible(model)
    transporter = {'FUMt2_2pp'};
    model.rev(ismember(model.rxns,transporter)) = 1;
    model.lb(ismember(model.rxns,transporter)) = -1000;
end

function model = makeIrreversible(model, rxns)
% make possible loop reactions irreversible
    if ~isempty(length(rxns))
        model.rev(ismember(model.rxns,rxns)) = 0;
        model.lb(ismember(model.rxns,rxns)) = 0;
    end
end