function runOptSwapB

    
    % options
    thKO = true;
    options.filename = 'g-coupled-import_all';
    
    
    %%%%%%%%%%%%%%%%%%%%%
    % setup
    setupCobraSolver;
    tw = setupTwitty; 
    
    % load model
    model = loadModelNamed('iAF1260b');
    thRxns = {'NADTRHD', 'THD2pp'};

    
    % load designs
    [numeric,txt,raw]=xlsread('g-coupled-import.xls','input');

    cols = raw(1,:);
    raw = raw(2:end,:);
    parameterCell = cell(0,5);
    
    % h = waitbar(0,'thRun');

    % output = cell(size(raw,2),8);
    % fullSolns = cell(size(raw,2),1);
    for i=1:size(raw,1)

        
        knockoutRxns = {};
        knockoutStr = '';
        for k=6:size(raw,2)
            if (~isnan(raw{i,k}))
                knockoutRxns = [knockoutRxns raw(i,k)];
                knockoutStr = [knockoutStr ' ' raw{i,k}];
            end
        end

        targetRxn = raw{i,4};
        exchangeRxn = raw{i,3};
        isAerobic = strcmp(raw(i,2),'aerobic');
        
        parameterCell(end+1,:) = {isAerobic, thKO, exchangeRxn, targetRxn, knockoutRxns};
  
    end
    

    optSwapB(model, thRxns, parameterCell, options)

    % save_to_base(1);
    tw.updateStatus(sprintf('@zakandrewking %s finished', options.filename)); 

end