%parameters
optFlag=0;
robustFlag=1;
knockNum=3;
maxWrobust=1e7;
objectiveInd = 150; %biomass
chemicalInd = 329; % ethanol
fixedGlu = 10;

%constants
maxWopt=1000;
if (robustFlag)
    maxW=maxWrobust;
else
    maxW=maxWopt;
end

results =  robustKnock(chemicalInd, objectiveInd,...
           knockNum, maxW,... 
           optFlag, robustFlag, fixedGlu);

 
 
 
 
  