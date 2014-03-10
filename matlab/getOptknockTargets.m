function [selectedRxns,sets,trimmedSets,geneList,single_KO_GRs,modelRed] = getOptknockTargets(model,ssExcludeList,nCarbonThr,shouldReduceModel)

%getOptknockTargets determines a reduced set of reactions which can be
%input into OptKnock for targets of deletion
%
% [selectedRxns,sets,trimmedSets,geneList,single_KO_GRs,modelRed] = getOptknockTargets(model,ssExcludeList,nCarbonThr)
%
% model            Structure containing all necessary variables to described a
%                  stoichiometric model
% ssExcludeList    List of subsystems that should NOT be included in the
%                  final list of reaction targets
% nCarbonThr       defines the min # of carbons that a metabolite---that is
%                  acted on in a reaction---can have in the final list of reactions
%
% selectedRxns     The selected targets for deletion in OptKnock
% single_KO_GRs    The growth rates (GRs) for each knock out (KO) examined
% sets             see 'identifyCorrelSets' function
% modelRed         see 'reduceModel' function
%
%
% Markus Herrgard and Adam Feist 2/7/07
% Updated, Zachary King 9/10/12

%%
% %1. reduce model
%[minflux, maxflux] = fastFVA(model, 0.1, 'max', 'cplex');
%modelRed = model;
%for i=1:size(model.rxns)
%    modelRed = changeRxnBounds(modelRed, model.rxns(i), minflux(i), 'l');
%    modelRed = changeRxnBounds(modelRed, model.rxns(i), maxflux(i), 'u');
%    if abs(minflux(i)) < 1e-9 && abs(maxflux(i)) < 1e-9    
%        modelRed = removeRxns(modelRed, modelRed.rxns(i));
        %modelRed = changeRxnBounds(modelRed, model.rxns(i), 0, 'b');
%    end
%end
%for i=1:size(model.rxns)
%    modelRed = changeRxnBounds(modelRed, model.rxns(i), minflux(i), 'l');
%    modelRed = changeRxnBounds(modelRed, model.rxns(i), maxflux(i), 'u');
%end

% run reduced model
if shouldReduceModel
    [modelRed,hasFlux,maxes,mins] = reduceModel(model,1e-9,false,false,true, ...
                                                true,true);
else
    modelRed = model;
end



%%
%2. remove the essential reactions (or those less than 5% of WT) from the list

%run simulation - get the WT GR
FBAsolution = optimizeCbModel(modelRed);
WT_max_GR = FBAsolution.f;
%save the original bounds
og_lb = modelRed.lb;
og_ub = modelRed.ub;
% find the exchange reactions
[selExc,selUpt] = findExcRxns(modelRed,1,1);

%initialize the desired reactions vector and KO GR vector
desired_rxn_set = selExc;


% Optional: Remove the reactions that are subjectively deemed not to be valid
% OptKnock targets
special_rxns = {...
%     'GLYCLTt2rpp', 'GLYCLTtex',... % OptKnock generates invalid results if these are KOed
    'GLCP2', 'GLCP', 'GLBRAN2', 'GLDBRAN2'...  % deal with glycogen
%     'FLDR', 'RNDR3b','RNDR3', 'RNDR4', 'RNDR4b', 'RNTR3c', 'RNTR4c', 'TRDR'...
    };  % contain different energy compounds than ATP, NAD(P)

% put both of the non-target reactions together in the string
non_target_rxns = {special_rxns{:}};


rxnID = findRxnIDs(modelRed,non_target_rxns);
sorted_rxnID = sort(rxnID); % gets rid of zeros in the ID vector
non_target_rxnID = sorted_rxnID(sorted_rxnID ~= 0); % gets rid of zeros in the ID vector

for b=1:length(non_target_rxnID)
    desired_rxn_set(non_target_rxnID(b),1) = true;
end

% computationally essential reactions
single_KO_GRs=zeros(length(modelRed.rxns),1);
for i=1:length(modelRed.rxns)
    if desired_rxn_set(i) == false % this skips all of the exchange and demand reactions
        %change the bounds
        modelRed = changeRxnBounds(modelRed,modelRed.rxns(i),0,'b');
        %run simulation
        FBAsolution = optimizeCbModel(modelRed);
        single_KO_GRs(i)=FBAsolution.f; % save the solution, note that EX_ and DM_ reactions are not tested.

        % mark the spot in the essential reactions vector if the growth rate is less than 5% of the WT GR
        if (FBAsolution.f/WT_max_GR) < 0.05
            desired_rxn_set(i) = true;  % almost essential reaction
        end
        %reset the model bounds
        modelRed = changeRxnBounds(modelRed,modelRed.rxns(i),og_lb(i),'l');
        modelRed = changeRxnBounds(modelRed,modelRed.rxns(i),og_ub(i),'u');
    end
end

%%
% generate 'selectedRxns'
s=1;
selectedRxns_v0={};
for r=1:length(modelRed.rxns)
    if desired_rxn_set(r) == false
        selectedRxns_v0(s,1)=modelRed.rxns(r);
        s=s+1;
    end
end

%%
%3. remove non-gene-associated reactions
u=1;
selectedRxns_v1={};
for t=1:length(selectedRxns_v0)
    if ~isequal(model.rules(find(strcmp(model.rxns,selectedRxns_v0(t,1)))),{''}) %looks to see if the gene association is empty
        selectedRxns_v1(u,1)=selectedRxns_v0(t,1);
        u=u+1;
    end
end

%%
%4. remove reactions from particular subsystem

v=1;
selectedRxns_v2={};
for w=1:length(selectedRxns_v1)
    if ~ismember(model.subSystems (find(strcmp(model.rxns,selectedRxns_v1(w,1) ))),ssExcludeList)
        selectedRxns_v2(v,1)=selectedRxns_v1(w,1);
        v=v+1;
    end
end

%%
% 5. remove reactions that act on molecules with > some number of carbons

% call the find carbon reactions
[hiCarbonRxns,zeroCarbonRxns,nCarbon] = findCarbonRxns(model,nCarbonThr); %used the whole model, not reduced
%join the two lists
hiAndZeroCarbonRxns = [hiCarbonRxns; zeroCarbonRxns];
%find which reactions in the list act on high carbon compounds
ref_CarbRxns = ismember(selectedRxns_v2,hiAndZeroCarbonRxns);
% generate the list
selectedRxns_v3 = selectedRxns_v2(~ref_CarbRxns);

%%
% 6. remove any reactions that are assigned to the 'SPONTANEOUS' gene: 's0001',
% fake gene product (it was generated for E coli for spontaneous and diffusion
% reactions)
[isInModel,geneInd] = ismember('s0001',model.genes);

if(isInModel)
    %find the indecies of reactions associated to the spontaneous genes
    rxnIndSpon = find(any(model.rxnGeneMat(:,geneInd),2));
else
    rxnIndSpon = [];
end
%find which reactions in the list are associated to SPONTANEOUS
ref_SponRxns = ismember(selectedRxns_v3,model.rxns(rxnIndSpon));
selectedRxns_v4 = selectedRxns_v3(~ref_SponRxns);

%%
% Optional: Remove the reactions that are subjectively deemed not to be valid
% OptKnock targets

%%%% need to write this if desired %%%%%

%%
%7. determine coupled reactions and only look at one, report the correlated
%sets

% %sampling with the existing reduced model
% warmupPts = createHRWarmup(modelRed,2000,true);
% ACHRSampler(modelRed,warmupPts,'Redmodel_samples_ROCT',1,100,10); % a relatively small number of points
% samples = loadSamples('Redmodel_samples_ROCT',1,100,0);

% samples based on the nullspeace
samples = null(full(modelRed.S));

%find the correlated sets
[sets,setNumber,setSize] = identifyCorrelSets(modelRed,samples);

% make 'trimmedSets_v0' of reactions that only contain reactions that are
% still considered targets up to this point
trimmedSets_v0={};
z=0;
for j = 1:length(sets) % each set
    y=0;
    for k = 1:length(sets{j,1}.set) % each element
        if ismember(sets{j,1}.names(k,1),selectedRxns_v4) == true; % is the reaction in the selected reastion? yes
            y=y+1;
            if y == 1
                z=z+1;
            end
            trimmedSets_v0{z,1}.set(y,1)=sets{j,1}.set(k,1);
            trimmedSets_v0{z,1}.names(y,1)=sets{j,1}.names(k,1);
        end
    end
end

%make a logical vector to select which reactions to keep as targets
desired_rxn_set2 = false(length(modelRed.rxns),1);
rxnID2 = findRxnIDs(modelRed,selectedRxns_v4);
for d=1:length(rxnID2)
    desired_rxn_set2(rxnID2(d),1) = true;
end

g=0;
trimmedSets={};
% use the set data to filter out more reactions
for p = 1:length(trimmedSets_v0) % each set
    if length(trimmedSets_v0{p,1}.set) > 1  %only if there is more than one element in the set
        g=g+1;
        trimmedSets{g,1}.set = trimmedSets_v0{p,1}.set;
        trimmedSets{g,1}.names = trimmedSets_v0{p,1}.names;
        for q=2:1:length(trimmedSets_v0{p,1}.set) %all of the elements in the set (except for the first one)
            desired_rxn_set2(trimmedSets_v0{p,1}.set(q,1)) = false; % only the first reaction in the set is still considered
        end
    end
end

%code to see how many reactions were removed by utilizing the co-sets
% count = 0;
% for i=1:length(trimmedSets)
% length(trimmedSets{i,1}.names)
% count = count + length(trimmedSets{i,1}.names);
% end

% generate 'selectedRxns' again
s=1;
selectedRxns={};
for r=1:length(modelRed.rxns)
    if desired_rxn_set2(r) == true
        selectedRxns(s,1)=modelRed.rxns(r);
        s=s+1;
    end
end
%%
% make a gene list of the genes that correspond to the selected reactions
geneList = false(size(model.genes));

for i = 1:length(selectedRxns)
    rxnindex = findRxnIDs(model, selectedRxns(i));
    geneList = geneList | model.rxnGeneMat(rxnindex,:)';
end

geneList = model.genes(geneList);