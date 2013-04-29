function [model, isSpecial] = setupModelForTarget(model, targetRxn)
    isSpecial = true;
    % deal with special cases
    if strcmp(targetRxn, 'EX_h2(e)') % turn on FHL for hydrogen production
        model = changeRxnBounds(model, 'FHL', 1000, 'u');
    % elseif strcmp(targetRxn, 'EX_ala-L(e)')  % exporting valine does not matter
    %     model = changeRxnBounds(model, 'EX_val-L(e)', 0, 'u');
    else
        isSpecial = false;
    end
end
