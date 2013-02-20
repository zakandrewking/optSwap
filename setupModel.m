function [model, biomassRxn] = ...
        setupModel(modelName,substrate,aerobicStr,transhydrogenaseKoStr,POR5Str)
% setupModel
% 
% INPUTS
%   modelName
%   substrate
%   aerobicStr - 'aerobic' or 'anaerobic'
%   transhydrogenaseKoStr - 'THKO' or 'noTHKO'
% 
% OUTPUTS
%   model
%   biomassRxn
%
% Zachary King 9/12/12


    
    aerobicStr = lower(aerobicStr);
    if strcmp(aerobicStr, 'aerobic')
        isAerobic = true;
    elseif strcmp(aerobicStr, 'anaerobic')
        isAerobic = false;
    else
        fprintf('misspelled %s\n', aerobicStr)
    end
    if nargin < 4
        transhydrogenaseKnockout = false;
    else
        transhydrogenaseKoStr = lower(transhydrogenaseKoStr);
        if strcmp(transhydrogenaseKoStr, 'thko')
            transhydrogenaseKnockout = true;
        elseif strcmp(transhydrogenaseKoStr, 'nothko')
            transhydrogenaseKnockout = false;
        else
            fprintf('misspelled %s\n', transhydrogenaseKoStr)
        end 
    end 
    if nargin < 5
        POR5Str = 'por5_irrev';
    else
        POR5Str = lower(POR5Str);
    end
    
    
    % load the model
    model = loadModelNamed(modelName);
    biomassRxn = model.rxns(model.c~=0);
    
    % set ATPM
    % model = changeRxnBounds(model, 'ATPM', 8.37, 'b');
    
    % prepare the model
    if transhydrogenaseKnockout
        % knock out th's
        thRxns = {'NADTRHD', 'THD2pp'};
        model = changeRxnBounds(model, thRxns, 0, 'b');
    end
    model = changeRxnBounds(model, 'EX_glc(e)', 0, 'l');
    if ~strcmp(substrate,'none')
        model = changeRxnBounds(model, substrate, -20, 'l');
    end
    if isAerobic
        model = changeRxnBounds(model, 'EX_o2(e)', -20, 'l');
    else
        model = changeRxnBounds(model, 'EX_o2(e)', 0, 'l'); 
        % BAD REACTIONS, DOWN BOY
        if strcmp(modelName,'iND750')
            necessary_ex = {'EX_ergst(e)', 'EX_zymst(e)', 'EX_hdcea(e)', ...
                            'EX_ocdca(e)', 'EX_ocdcea(e)', ...
                            'EX_ocdcya(e)'};
            % 'EX_k(e)','EX_na1(e)', 'EX_co2(e)'};
            model = changeRxnBounds(model, necessary_ex, -1000, 'l');
            model = changeRxnBounds(model, necessary_ex, 1000, 'u'); 
        else
            model = changeRxnBounds(model, {'CAT';'SPODM';'SPODMpp'}, [0;0;0], 'b');
        end
        % can make oxygen...turn off for anaerobic simulations 
    end
    
    % Make POR5 irreversible
    if ~strcmp(modelName,'iND750')
        if strcmp(POR5Str,'por5_irrev') 
            model = changeRxnBounds(model, 'POR5', 0, 'l');
            model.rev(ismember(model.rxns, 'POR5')) = 0;
        elseif strcmp(POR5Str,'por5_rev')
        else
            error(sprintf('misspelled %s', POR5Str));
        end
    end
    
end