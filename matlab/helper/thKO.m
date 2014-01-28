function [soln,thFlux,thFraction,nadphFlux] = ...
        thKO(model,knockoutRxns,targetRxn,biomassRxn,uptakeRxns,isAerobic,noEnv)

    % THKO
    %
    % [soln,thFlux,thFraction,nadphFlux]
    %
    % runs FBA of model and determines flux through transhydrogenase reactions thRxns
    %
    % INPUTS
    % model
    % knockoutRxns
    % targetRxn
    % biomassRxn
    %
    % OPTIONAL INPUTS
    % uptakeRxns - carbon source reactions
    % isAerobic
    % noEnv - don't make production envelope
    %
    % OUTPUTS
    %
    % sol - solutions to FBAs
    % thFlux - fluxes through transhydrogenase reactions
    % thFraction - fraction of nadph produced by th
    %


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PARAMETERS - set parameters here

    nreq = 4;

    if nargin < nreq + 1
        uptakeRxns = 'EX_glc(e)';
    end
    if nargin < nreq + 2
        isAerobic = false;
    end
    if nargin < nreq + 3
        noEnv = false;
    end

    uptakeMax = 20;
    O2Uptake = 20;
    thRxns = {'NADTRHD', 'THD2pp'};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Set up model bounds

    % set o2 bounds
    if isAerobic
        model = changeRxnBounds(model, 'EX_o2(e)', -O2Uptake, 'l');
    else
        model = changeRxnBounds(model, 'EX_o2(e)', 0, 'l');
    end

    % turn off glucose
    model = changeRxnBounds(model, 'EX_glc(e)', 0, 'l');

    % set carbon substrate bounds
    model = changeRxnBounds(model, uptakeRxns, -uptakeMax, 'l');

    % turn on FHL for H2 simulations
    if (strcmp(targetRxn,'EX_h2(e)'))
        model = changeRxnBounds(model, 'FHL', 1000, 'u');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Create production envelope

    if (~noEnv)

        display(knockoutRxns);

        figure(1);
        productionEnvelope(model,knockoutRxns, 'k',...
                           targetRxn,biomassRxn);

        xlabel('biomass output (hr^-^1)',...
               'FontSize', 14);
        ylabel([targetRxn ' production (mmol gDW^-^1 hr^-^1)'],...
               'FontSize', 14);

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% FBA and find flux through THs

    model = changeRxnBounds(model, knockoutRxns, 0, 'b');
    % model = changeObjective(model, targetRxn);
    %        dbstop if (strcmp(uptakeRxns,'EX_xyl_D(e)'))
    soln = optimizeCbModel(model);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Find flux through th's

    if (isempty(soln) || isempty(soln.x))
        display('empty flux vector')
        nadphFlux = 0;
        thFlux = 0;
        thFraction = 0;
    else
        nadphFlux = fluxMakingNADPH(model,soln.x);
        thFlux = fluxMakingNADPH(model,soln.x,thRxns);
        thFraction = thFlux/nadphFlux;
    end

    if (strcmp(targetRxn,'EX_h2(e)'))% && isAnaerobic)
        FHLflux = soln.x(ismember(model.rxns,'FHL'))
    end

end



function sumFlux = fluxMakingNADPH(model,x,rxns)

    theMet = 'nadph[c]';

    if nargin < 3
        rxns = model.rxns;
    end

    % reverse fluxes when nadph is a reactant
    metInd = find(ismember(model.mets,theMet));
    nadphFluxes = x.*model.S(metInd,:)';
    selectedFluxes = full(nadphFluxes(ismember(model.rxns,rxns)));

    % return total positive flux (ignore nadph consuming reactions)
    sumFlux = sum(selectedFluxes(selectedFluxes>0));

end