function runThKO()


    % cd ~/Dropbox/lab/OptSwap/code
    % fid = fopen('/Users/zaking/Dropbox/MATLABoutput/log.txt','at');
    
    setupCobraSolver;
    
    %%% load model
    model = loadModelNamed('iAF1260b');

    
    % load designs
    [numeric,txt,raw]=xlsread('g-coupled-import.xls','input');

    cols = raw(1,:);
    raw = raw(2:end,:);

    h = waitbar(0,'thRun');

    output = cell(size(raw,2),8);
    fullSolns = cell(size(raw,2),1);
    for i=1:size(raw,1)
        display(i)

        knockoutRxns = {};
        knockoutStr = '';
        for k=6:15
            if (~isnan(raw{i,k}))
                knockoutRxns = [knockoutRxns raw(i,k)];
                knockoutStr = [knockoutStr ' ' raw{i,k}];
            end
        end

        targetRxn = raw{i,4};
        exchangeRxn = raw{i,3};
        isAnaerobic = strcmp(raw(i,2),'anaerobic');
        
        switch isAnaerobic
          case true
            anString = 'anaerobic';
          case false
            anString = 'aerobic';
        end
        
        [sol,thFlux,thFraction,nadphFlux] = ...
            thKO(model,knockoutRxns,targetRxn,...
                 exchangeRxn,isAnaerobic,true);
        
       
        output(i,1:8) = {nadphFlux,thFlux,thFraction,sol.f,...
                         anString, exchangeRxn, targetRxn,...
                         knockoutStr};
        
        % fprintf(fid, [num2str(i) ' \n']);
        fullSolns{i} = sol;
        % gr = sol.f

        waitbar(i/size(raw,1),h);
    end
    
    cell2csv('out.csv',output)

    % fclose(fid);
    
    close(h)

    %save variables to workspace
    save_to_base(1);
end