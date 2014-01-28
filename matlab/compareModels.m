function output = compareModels(modelCell, opt)

    if nargin < 1, error('Not enough arguments'); end
    if ~iscell(modelCell), modelCell = {modelCell}; end
    if ~exist('opt','var')
        opt = struct();
    end
    if ~isfield(opt,'targetRxn'), opt.targetRxn = 'EX_etoh(e)'; end
    if ~isfield(opt,'substrateRxn'), opt.substrateRxn = 'EX_glc(e)'; end
    if ~isfield(opt,'isAerobic'), opt.isAerobic = false; end
    if ~isfield(opt,'O2Uptake'), opt.O2Uptake = 20; end
    if ~isfield(opt,'maxSubsUptake'),  opt.maxSubsUptake = 20; end
    if ~isfield(opt,'mapOutputFormat'), opt.mapOutputFormat = 'svg'; end
    if ~isfield(opt,'pairwiseCompare'), opt.pairwiseCompare = false; end
    if ~isfield(opt,'differenceThreshold'), opt.differenceThreshold =  1e-5;
    if ~isfield(opt,'showEnvelope'), opt.showEnvelope = false; end

    changeCbMapOutput(opt.mapOutputFormat);
    if isunix
        mapPath = '~/lab/iAF1260-map.txt';
    else
        mapPath = 'C:\Documents and Settings\z1king\My Documents\Dropbox\lab\iAF1260-map.txt';
    end
    
    map = readCbMap(mapPath);

    resultCell = {};
    for i = 1:length(modelCell)
        disp(sprintf('model %i', i));

        model = modelCell{i};

        % setup models
        if opt.isAerobic
            model = changeRxnBounds(model, 'EX_o2(e)',...
                                    -opt.O2Uptake, 'l');
        else
            model = changeRxnBounds(model, 'EX_o2(e)', 0, 'l');
        end

        model = changeRxnBounds(model, 'EX_glc(e)', 0, 'l');
        model = changeRxnBounds(model, opt.substrateRxn, ...
                                -opt.maxSubsUptake, 'l');

        result = optimizeCbModel(model);
        if isempty(result.x)
            warning('could not optimize');
        else
            targetProduction = result.x(ismember(model.rxns, opt.targetRxn))
            % disp(sprintf('target production: %g', targetProduction));

            if opt.showEnvelope
                figure(i)
                productionEnvelope(model, {}, '-k', opt.targetRxn, ...
                                   model.rxns(model.c~=0));
            end
            
            resultCell(end+1) = {result};
        end

    end

    difOutput = {};
    if opt.pairwiseCompare
        if size(resultCell{2}.x) ~= size(resultCell{1}.x)
            error('model rxn sets do not match')
        end

        dif = resultCell{2}.x - resultCell{1}.x;
        
        options.fileName = ['model.svg'];
        drawFlux(map, model, dif, options);
        
        th = opt.differenceThreshold;
        difCell = [num2cell(dif(dif<-th | dif>th)), model.rxns(dif<-th | dif>th), ...
            model.subSystems(dif<-th | dif>th)];
        difOutput{end+1} = sortcell(difCell, 1);
        disp(difOutput{1})
        output = difOutput{1};
    else
        output = {};
    end

end