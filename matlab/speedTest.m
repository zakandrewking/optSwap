function speedTest()
    logFile = 'database.csv';
    global fileId
    fileId = fopen(logFile, 'a');
    
    display('setting up model')
    % model = setupModel('iJO', 'EX_glc(e)', false, true);
    opt.targetRxn = 'EX_for(e)';
    % model = setupModelForTarget(model, opt.targetRxn);
    experiment = 'speedTest 10KO';
    run = {
           'optKnock-MixedSel-red', 'robustKnock-MixedSel-red',...
           'CbOptKnock-CbSel-red'...
          };
           % 'optKnock-CbSel','robustKnock-CbSel',...
           % 'optKnock-MixedSel', 'robustKnock-MixedSel',...
           % 'optKnock-RobustSel', 'robustKnock-RobustSel',...
           % 'optKnock-CbSel-red','robustKnock-CbSel-red',...
           % 'optKnock-RobustSel-red', 'robustKnock-RobustSel-red',...

    knockTypes = [0,1,0,1,0,1];
     
    opt.knockoutNum = 10;
    
    % Choose reactions for knocking
    % display('selecting reactions');
    % ssExcludeList = {'Cell Envelope Biosynthesis',...
    %                  'Exchange',...
    %                  'Inorganic Ion Transport and Metabolism',...
    %                  'Lipopolysaccharide Biosynthesis / Recycling',...
    %                  'Murein Biosynthesis',...
    %                  'Murein Recycling',...
    %                  'Transport, Inner Membrane',...
    %                  'Transport, Outer Membrane',...
    %                  'Transport, Outer Membrane Porin',...
    %                  'tRNA Charging',...
    %                 };
    % nCarbonThr = 10;
    % [selectedRxns,~,~,~,~,~] = ...
    %     getOptknockTargets(model,ssExcludeList,nCarbonThr,false);
    
    % rTic = tic;
    % [selectedRxnsRed,~,~,~,~,modelRed] = ...
    %     getOptknockTargets(model,ssExcludeList,nCarbonThr,true);
    % save('reducedModel-iJO-for.mat', 'selectedRxnsRed','modelRed');
    load('reducedModel-iJO-for.mat')
    tRed = 0;
    % tRed = toc(rTic)/60;
    % display(sprintf('%.1f min to reduce model & select reactions,', tRed));
    
    % outputSelectedRxns = cell(1,1);
    for i=1:2
        lTic = tic;
        display(sprintf('set %d: %s', i, run{i}));
        % if i<=6
        %     modelT = model;
        %     selectedRxnsT = selectedRxns;
        % else
            modelT = modelRed;
            selectedRxnsT = selectedRxnsRed;
        % end
        % switch i
        %   case {1,2}
        %     opt.knockableRxns = selectedRxnsT;
        %     opt.notKnockableRxns = modelT.rxns(~ismember(modelT.rxns, ...
        %                                                 selectedRxnsT));
        %   case {3,4}
            opt.knockableRxns = {};
            opt.notKnockableRxns = modelT.rxns(~ismember(modelT.rxns, ...
                                                        selectedRxnsT));
        %   case {5,6}
        %     opt.knockableRxns = {};
        %     opt.notKnockableRxns = {};
        % end
        
        opt.knockType = knockTypes(i);
        
        results = optSwapD(modelT, opt);
        % outputSelectedRxns{i} = modelT.rxns(results.yInd);
        
        
        t = toc(lTic);
        
        myPrint('OptSwap,', []);
        myPrint('%s,', experiment);
        myPrint('%s,', datestr(now, 'mm/dd/yyyy'));
        myPrint('%s,', datestr(now, 'HH:MM:SS'));
        myPrint('%s,', run{i});
        myPrint('optSwapD.m,', []);
        myPrint('%d,', opt.knockoutNum);
        myPrint('0,', []);
        myPrint('%s,', opt.targetRxn);
        for j=1:4, myPrint(',', []); end
        if ~isempty(results.knockoutRxns)
            printReactions(results.knockoutRxns);
        else
            display('no knockouts');
            myPrint(',', []);
        end
        if ~isempty(results.knockoutDhs)
            printReactions(results.knockoutDhs);
        else
            display('no swaps');
            myPrint(',', []);
        end
        myPrint('%.1f,', t/60);
        myPrint(',', []);
        myPrint(',', []);
        myPrint('%.1f min to reduce model,', tRed);
        myPrint('\n', []);
    end
    
    % save('outputSelectedRxns.mat', 'outputSelectedRxns');
    
    for i=3:3
x        % switch i
        %     case 13
        %       modelT = model;
        %       selectedRxnsT = selectedRxns;
        %     case 14
              modelT = modelRed;
              selectedRxnsT = selectedRxnsRed;
        % end
        lTic = tic;
        options.numDel = opt.knockoutNum;
        options.targetRxn = opt.targetRxn;
        optKnockSol = OptKnock(modelT,selectedRxnsT,options);
        t = toc(lTic);
        myPrint('OptSwap,', []);
        myPrint('%s,', experiment);
        myPrint('%s,', datestr(now, 'mm/dd/yyyy'));
        myPrint('%s,', datestr(now, 'HH:MM:SS'));
        myPrint('%s,', run{i});
        myPrint('OptKnock.m,', []);
        myPrint('%d,', opt.knockoutNum);
        myPrint('0,', []);
        myPrint('%s,', opt.targetRxn);
        for j=1:4, myPrint(',', []); end
        if ~isempty(optKnockSol.rxnList)
            printReactions(optKnockSol.rxnList);
        else
            display('no knockouts');
            myPrint(',', []);
        end
        myPrint(',', []);
        myPrint('%.1f,', t/60);
        myPrint(',', []);
        myPrint(',', []);
        myPrint('%.1f min to reduce model,', tRed);
        myPrint('\n', []);
    end
    
    % tw = setupTwitty; 
    % tw.updateStatus(sprintf('@zakandrewking --> speedTest finished'));
end

function printReactions(reactions)
    for j=1:length(reactions)
        myPrint('''%s'';', reactions{j});
    end
    myPrint(',', []);
end

function myPrint(string, val)
    global fileId
    display(sprintf(string, val));
    fprintf(fileId, string, val);
end
