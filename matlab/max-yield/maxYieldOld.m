function maxYieldOld
% maxYield
% calculate maximum yield of target reactions at minimum biomass production
% with dehydrogenase swaps
% 
% Zachary King 9/20/2012
    
    cleaner = onCleanup(@() cleanup);
    global run status
    status = 'starting';
    run = 'max-yield_all-forward-transporters';
    
    logFile = sprintf('%s-%s.csv', run, datestr(now, 'yy-mm-dd_HH_MM_SS'));
    fileId = fopen(logFile, 'a');
    fprintf('%s\n', logFile);

    % load model
    model = loadModelNamed('iJO');
    [model, biomassRxn] = setupModel('iJO','EX_glc(e)', 'anaerobic', ...
                                           'noTHKO');
    model.c(:) = 0;
    
    % get all exchange reactions
    targetRxns = model.rxns(findExcRxns(model));
    % turn on special cases
    model = makeFumTransporterReversible(model);
    % remove exchanges that don't exit cell
    noExitRxns = model.rxns(model.ub <= 0);
    targetRxns(ismember(targetRxns,noExitRxns)) = [];
    % set parameters 
    dhRxns = dhRxnList(31);
    growthMins = [0.1];
    isAerobic = false;
    substrateList = {'EX_glc(e)', 'EX_xyl-D(e)'};

    % testing
    % targetRxns = {'EX_for(e)'};
    
    % go!
    out = zeros(3, length(substrateList), length(targetRxns), length(growthMins));
    fprintf(fileId, 'swap,substrate,target,growth min,max yield\n'); 
    for m = 1:3 
        model_K = model;
        switch m
          case 1
            swapStr = 'wildtype';
          case 2
            swapStr = 'THKO';
            % knock out transhydrogenases
            thRxns = {'NADTRHD', 'THD2pp'};
            model_K = changeRxnBounds(model_K, thRxns, 0, 'b');
          case 3
            swapStr = 'DH-swap';
            thRxns = {'NADTRHD', 'THD2pp'};
            model_K = changeRxnBounds(model_K, thRxns, 0, 'b');
            [model_K, newNames] = modelSwap(model_K, dhRxns, true);
        end

        for i=1:length(substrateList) 
            substrate = substrateList{i}; 
            model_KS = model_K;
            model_KS = changeRxnBounds(model_KS, 'EX_glc(e)', 0, 'l');
            model_KS = changeRxnBounds(model_KS, substrate, -20, 'l');

            for j = 1:length(targetRxns) 
                targetRxn = targetRxns{j};
                model_KST = model_KS;
                model_KST = setupModelForTarget(model_KST, targetRxn);
                model_KST.c = zeros(size(model_KST.c));
                model_KST.c(ismember(model_KST.rxns, targetRxn)) = 1;
                
                for k = 1:length(growthMins) 
                    growthMin = growthMins(k);
                    model_KSTG = model_KST; 
                    model_KSTG.lb(ismember(model_KSTG.rxns, biomassRxn)) = growthMin;

                    disp('running optimization');
                    status = sprintf('%d, %s, %s, min: %.1f', m, substrate, targetRxn, growthMin); 
                    sol = optimizeCbModel(model_KSTG); 
                    out(m,i,j,k) = sol.f; 
                    fprintf(fileId, '%s,%s,%s,%.1f,', swapStr, substrate, ...
                            targetRxn, growthMin);
                    fprintf(fileId, '%.2f,\n', sol.f);
                    % check FUM flux
                    % if ~isempty(sol.x)
                    %     fprintf('Fum flux: %.2f\n', ...
                    %             sol.x(ismember(model.rxns,'FUMt2_2pp')));
                    % end
                end
            end
        end
    end
    status = 'finished'
    fclose(fileId);
    save('raw.mat', 'out');
end

function model = makeFumTransporterReversible(model)
    transporter = {'FUMt2_2pp'};
    model.rev(ismember(model.rxns,transporter)) = 1;
    model.lb(ismember(model.rxns,transporter)) = -1000;
end





                    % if ~isempty(sol.x)
                        % fluxes{i,j,k,m} = sol.x(ismember(model_KSTG.rxns,dhRxns));
                        % if m==3
                        %     swapFluxes{i,j,k,m} = sol.x(ismember(model_KSTG.rxns, ...
                        %                                          newNames));
                        % end
                    % else
                    %     fluxes{i,j,k,m} = [];
                    % end 

    
    % global fileId dhRxns newNames targetRxns growthMins fluxes swapFluxes
    % fluxes = cell(size(out)); swapFluxes = fluxes;
    % fprintf(fileId, '\n\n\nFLUXES\n\n');
      
    % for i=1:length(substrateList)
    %     printSubstrate(substrateList{i},i);
    % end




% function printSubstrate(substrate,i)
%     global fileId targetRxns
%     fprintf(fileId, '%s\n\n', substrate);

%     for j=1:length(targetRxns)
%         printTarget(targetRxns{j},i,j);
%     end

%     fprintf(fileId, '\n');
% end

% function printTarget(target,i,j)
%     global fileId dhRxns newNames
%     fprintf(fileId, 'target: %s\n', target);
%     fprintf(fileId, ',0.0,,,0.1,,,0.2,,,0.3\n,');
%     for x=1:4
%         fprintf(fileId, 'noThKO,THKO,THKO+swaps,');
%     end
%     fprintf(fileId, '\n');
%     for n=1:length(dhRxns)
%         printDhRxn(dhRxns{n},i,j,n,false);
%     end
%     for n=1:length(newNames)
%         printDhRxn(newNames{n},i,j,n,true);
%     end
% end

% function printDhRxn(dhRxn,i,j,n,swap)
%     global fileId growthMins fluxes swapFluxes
%     fprintf(fileId, '%s,', dhRxn);
%     for k=1:length(growthMins)
%         for m=1:3
%             if swap && (m < 3)
%                 fprintf(fileId, ',');
%                 continue;
%             elseif swap
%                 fl = swapFluxes{i,j,k,m};
%             else
%                 fl = fluxes{i,j,k,m};
%             end
%             if isempty(fl), fl = 0;
%             else fl = fl(n); end
%             fprintf(fileId, '%.2f,', fl);
%         end
%     end
%     fprintf(fileId, '\n');
% end
