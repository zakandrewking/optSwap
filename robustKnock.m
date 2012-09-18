function results=robustKnock(chemicalInd, objectiveInd, knockoutNum, maxW, optKnockFlag, robustKnockFlag, fixedGlu)
%chemicalInd - chemical to produce
%objectiveInd - organizms objective
%knockoutNum - number of knockouts
%maxW - maximal value of dual variables (higher number will be more
%accurate but takes more calculation time)
%optKnockFlag - indicates if optKnock will be calculated
%robustKnockFlag - indicates if optKnock will be calculated
%fixedGlu - size of fixed glucose reaction

%parameters
findMaxWFlag=0;
P=1;
coupledFlag = 1;

K=knockoutNum;
consModel=createModel(chemicalInd, objectiveInd, fixedGlu);
[model,yInd, notYInd, m,n,notY, coupled, coupledYInd, coupledNotYInd]=...
    createRobustOneWayFluxesModel(consModel, chemicalInd, coupledFlag);

ySize=size(yInd,1);

%part 2
%max C'v
%s.t
%[A,Ay]*[v;y]<=B
coupledYsize=size(coupledYInd,1);
I=eye(m);
A=[ model.S;...
   -model.S;...
    I(notYInd,:);...
    I(coupledNotYInd,:);...
   -I(notYInd,:);...
   -I(coupledNotYInd,:);...
    I(yInd,:);...
    I(coupledYInd,:);...
   -I(yInd,:);...
   -I(coupledYInd,:)];

[aSizeRow, vSize]=size(A);
knocables=zeros(m,1);
knocables(yInd)=1;

Ay1=diag(knocables);
Ay1(coupled(:,2), :)=Ay1(coupled(:,1), :);
Ay1=Ay1*diag(model.ub);
for j=1:length(coupled)
    if (model.ub(coupled(j,1)) ~=0)
        Ay1(coupled(j,2), coupled(j,1))=Ay1(coupled(j,2), coupled(j,1)).*(model.ub(coupled(j,2))./model.ub(coupled(j,1)));
    end
end

Ay2=diag(knocables);
Ay2(coupled(:,2), :)=Ay2(coupled(:,1), :);
Ay2=Ay2*diag(model.lb);
for j=1:length(coupled)
    if (model.lb(coupled(j,1)) ~=0)
        Ay2(coupled(j,2), coupled(j,1))=Ay2(coupled(j,2), coupled(j,1)).*(model.lb(coupled(j,2))./model.lb(coupled(j,1)));
    end
end

z1=find(Ay1);
z2=find(Ay2);
zSize=size([z1;z2],1);

Ay=[zeros(2*n+2*(vSize-ySize-coupledYsize),ySize);
    -Ay1(yInd,yInd);
    -Ay1(coupledYInd,yInd);
    Ay2(yInd,yInd);
    Ay2(coupledYInd,yInd);];  %flux boundry constraints

%so: [A,Ay]x<=B;
B=[zeros(2*n,1);
   model.ub(notYInd);
   model.ub(coupledNotYInd);
   -model.lb(notYInd);
   -model.lb(coupledNotYInd);
   zeros(2*(ySize+coupledYsize),1);    ];

C=model.organismObjective;

%needs to be minimum so:
%min -C*v
%s.t
%-[A,Ay]*[v;y]>=-B
[A_w, Ay_w ,B_w,C_w, lb_w, ub_w, wSize, wZs]=seperateTransposeJoin(-A, -Ay,-B,-C,ySize, 1,  m, maxW,findMaxWFlag, zSize);
%max C_w(w  z)'
%s.t
%[A_w  Ay_w]*(w z y)  <=  B_w
awSizeRow=size(A_w,1);

Ajoined=[
%C', P*C_w', sparse(1, ySize);                                          %dual problem objective function = primal problem objective function
    -C', -P*C_w', sparse(1, ySize);
    A, sparse(aSizeRow, wSize+zSize), Ay;               %stochiometric constraints + flux boundry constraints
    sparse(awSizeRow, vSize), A_w, Ay_w;                  %dual constraints
    zeros(1,vSize+wSize+zSize), -ones(1,ySize);
        ];

Bjoined=[
    0;
% 0;
    B;
    B_w;
    K-ySize;
        ];

Cjoined=[model.C_chemical; zeros(wSize+zSize,1); zeros(ySize,1)];

tmpLoptKnock=lb_w(1:wSize+zSize);
tmpHoptKnock=ub_w(1:wSize+zSize);
ysUpperBoundOptKnock=ones(ySize,1);

lbJoined=[model.lb;tmpLoptKnock; zeros(ySize, 1)];
ubJoined=[model.ub;tmpHoptKnock; ysUpperBoundOptKnock];

%max Cjoined*x'
%s.t
%Ajoined*x  <=  Bjoined

IntVars_optKnock=vSize+wSize+zSize+1:vSize+wSize+zSize+ySize;

%solve milp
%parameter for mip assign
x_min = []; x_max = []; f_Low = -1E7; % f_Low <= f_optimal must hold
f_opt = -141278;
nProblem = 7; % Use the same problem number as in mip_prob.m
fIP = []; % Do not use any prior knowledge
xIP = []; % Do not use any prior knowledge
setupFile = []; % Just define the Prob structure, not any permanent setup file
x_opt = []; % The optimal integer solution is not known
VarWeight = []; % No variable priorities, largest fractional part will be used
KNAPSACK = 0; % First run without the knapsack heuristic

if (optKnockFlag ==1)
    Prob_optKnock=mipAssign(-Cjoined, Ajoined, [], Bjoined, lbJoined, ubJoined, [], 'part 2.3 MILP', ...
                            setupFile, nProblem, ...
                            IntVars_optKnock, VarWeight, KNAPSACK, fIP, xIP, ...
                            f_Low, x_min, x_max, f_opt, x_opt);
    Prob_optKnock=setParams(Prob_optKnock);
    Result_optKnock=tomRun('cplex', Prob_optKnock, 1);

    YsOptKnockNotRound=Result_optKnock.x_k(IntVars_optKnock);
    YsOptKnock=round(YsOptKnockNotRound);
    knockedYind=find(1.-YsOptKnock);
    knockedVoptKnock=yInd(knockedYind);

    results.Result_optKnock.exitFlag=Result_optKnock.ExitFlag;
    results.Result_optKnock.organismObjectiveInd=model.organismObjectiveInd;
    results.Result_optKnock.chemicalInd=model.chemicalInd;
    results.Result_optKnock.objective=Result_optKnock.x_k(model.organismObjectiveInd);
    results.Result_optKnock.chemicel=Result_optKnock.x_k(model.chemicalInd);
    results.Result_optKnock.x_k=Result_optKnock.x_k;
    results.Result_optKnock.ySize=ySize;
    results.Result_optKnock.ySum=sum(YsOptKnockNotRound);
    results.Result_optKnock.yBoundry=ySize-K;
    results.Result_optKnock.knockedYind=knockedYind;
    results.Result_optKnock.knockedYvals=YsOptKnockNotRound(knockedYind);
    results.Result_optKnock.knockedV=knockedVoptKnock;
    results.Result_optKnock.knockedVvalues=Result_optKnock.x_k(yInd(results.Result_optKnock.knockedYind));
    results.Result_optKnock.y=YsOptKnockNotRound;

    if(Result_optKnock.ExitFlag ==0)        %a sucessful run
        results.Result_optKnock.knockCheck=knockCheck(model, knockedVoptKnock, coupledFlag, coupled, maxW);
    end
end

if (robustKnockFlag)
    %**************************************************************************
    %part 3

    %max min Cjoined*x'
    %s.t
    %Ajoined*x  <=  Bjoined

    %min Cjoined*x'
    %s.t
    %Ajoined*x  <=  Bjoined
    %equals to
    %min Cjoined*x'
    %s.t
    %-Ajoined*x  >=  -Bjoined

    A2=-[
    %C', P*C_w';
        -C', -P*C_w';
        A, sparse(aSizeRow, wSize+zSize);
        sparse(awSizeRow, vSize), A_w];

    Ay2=-[
    %zeros(1, ySize);
        zeros(1, ySize);
        Ay;
        Ay_w];

    C2=Cjoined(1:vSize+wSize+zSize,:);
    B2=-[
    % 0;
        0;
        B;
        B_w;
        ];

    z3=find(Ay2);
    zSizeOptKnock2=size(z3,1);

    [A2_w, Ay2_w ,B2_w,C2_w, lb2_w, ub2_w, uSize, uZs]=...
        seperateTransposeJoin(A2, Ay2,B2,C2 ,ySize, 1,  vSize+wSize+zSize,...
                              maxW,findMaxWFlag, zSizeOptKnock2);

    %max C2_w*x'
    %s.t
    %A2_w*x+Ay2_w*y  <=  B2_w

    %add u1, u2 variables so y will be feasible. add the constraints:
    % su=0 and umin*y<u<umax*y
    [A2_wRow, A2_wCol]=size(A2_w);
    [ARow, ACol]=size(A);

    A3=[
        A2_w, sparse(A2_wRow,ACol), Ay2_w;
           %dual constraints
        zeros(1,uSize+zSizeOptKnock2+vSize),  -ones(1,ySize);      
           %y sum constraints
        sparse(ARow,uSize+zSizeOptKnock2), A, Ay                             
           %feasibility conatraint
       ];

    B3=[
        B2_w;
        K-ySize;
        B
       ];

    C3=[C2_w;
        zeros(ACol,1);
        zeros(ySize,1)
       ];

    tmpL=lb2_w(1:uSize+zSizeOptKnock2);
    tmpH=ub2_w(1:uSize+zSizeOptKnock2);
    ysUpperBound=ones(ySize,1);
    lb3=[tmpL; model.lb; zeros(ySize,1)];
    ub3=[tmpH; model.ub; ysUpperBound];
    intVars=A2_wCol+ACol+1:A2_wCol+ACol+ySize;

    %solving mip
    Prob_OptKnock2=mipAssign(-C3, A3, [], B3, lb3, ub3, [], 'part 3 MILP', ...
                             setupFile, nProblem, ...
                             intVars, VarWeight, KNAPSACK, fIP, xIP, ...
                             f_Low, x_min, x_max, f_opt, x_opt);

    Prob_OptKnock2=setParams(Prob_OptKnock2);

    Result_OptKnock2=tomRun('cplex', Prob_OptKnock2, 2);

    YsNotRound=Result_OptKnock2.x_k(intVars);
    Ys=round(YsNotRound);
    knockedYindNotRound=find(1.-YsNotRound);
    knockedYindRound=find(1.-Ys);
    knockedVindRobustKnock=yInd(knockedYindRound);

    results.Result_OptKnock2.exitFlag=Result_OptKnock2.ExitFlag;
    results.Result_OptKnock2.organismObjectiveInd=model.organismObjectiveInd;
    results.Result_OptKnock2.chemicalInd=model.chemicalInd;
    results.Result_OptKnock2.chemical=-Result_OptKnock2.f_k;
    results.Result_OptKnock2.exitFlag=Result_OptKnock2.ExitFlag;
    results.Result_OptKnock2.ySize=ySize;
    results.Result_OptKnock2.ySum=sum(Result_OptKnock2.x_k(intVars));
    results.Result_OptKnock2.yLowerBoundry=ySize-K;
    results.Result_OptKnock2.knockedYindRound=knockedYindRound;
    results.Result_OptKnock2.knockedYvalsRound=YsNotRound(knockedYindRound);
    results.Result_OptKnock2.y=Ys;
    results.Result_OptKnock2.knockedVind=knockedVindRobustKnock;

    results.Result_OptKnock2.knockCheck = ...
        knockCheck(model, knockedVindRobustKnock, coupledFlag, coupled, maxW);

end

% fileName2=sprintf('res_%s.mat', datestr(clock, 'yyyymmddTHHMMSS'));
% fl = sprintf('%s//%s', pth, fileName2);
%save (fl);



