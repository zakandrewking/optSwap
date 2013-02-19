function results = optSwap(model, opt)
% Adapted from RobustKnock.m
% By Zachary King 8/13/2012
%
% INPUTS
% model
%
% OPTIONAL
% opt - struct with any of the following options:
%   knockType - 2 run OptSwap; 1 run robustKnock
%   targetRxn - chemical to produce
%   swapNum - number of swaps (-1 no limit)
%   knockoutNum - number of knockouts (-1 no limit)
%   interventionNum - maximum number of interventions (-1 no limit)
%   knockableRxns - reactions that should definitely be considered for knockout
%   notKnockableRxns - reactions that should be excluded from knockout
%   useCobraSolver - 1 use any cobra solver; 0 use tomlab cplex solver
%   maxW - maximal value of dual variables (higher number will be more
%          accurate but takes more calculation time)
%   biomassRxn - biomass objective function in the model
%   dhRxns - oxidoreductase reactions that can be swapped
%   maxTime - time limit in minutes
%   allowDehydrogenaseKnockout - use less than or equal constraint to allow
%                                dehydrogenase reaction knockouts
%
% OUTPUTS
% results

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up parameters
    if nargin < 1, error('Not enough arguments'); end

    % set default parameters
    if ~exist('opt','var')
        opt = struct();
    end
    if ~isfield(opt,'knockType'), opt.knockType = 0; end
    if ~isfield(opt,'targetRxn'), opt.targetRxn = 'EX_for(e)'; end
    if ~isfield(opt,'swapNum'), opt.swapNum = 0; end
    if ~isfield(opt,'knockoutNum'), opt.knockoutNum = 3; end
    if ~isfield(opt,'interventionNum'), opt.interventionNum = -1; end
    if ~isfield(opt,'knockableRxns'), opt.knockableRxns = {}; end
    if ~isfield(opt,'notKnockableRxns'), opt.notKnockableRxns = {}; end
    if ~isfield(opt,'useCobraSolver'), opt.useCobraSolver = 0; end
    if ~isfield(opt,'maxW'),
        switch (opt.knockType)
          case 2, opt.maxW = 1e7;
          case 1, opt.maxW = 1e7;
          case 0, opt.maxW = 1000;
        end
    end 
    if ~isfield(opt,'biomassRxn')
        opt.biomassRxn = model.rxns(model.c~=0);
    end
    if ~isfield(opt,'dhRxns'), opt.dhRxns = {}; end
    if ~isfield(opt, 'maxTime'), opt.maxTime = []; end
    if ~isfield(opt, 'printIntermediateSolutions')
        opt.printIntermediateSolutions = false;
    end
    if ~isfield(opt, 'intermediateSolutionsFile')
        opt.intermediateSolutionsFile = [];
    end
    if ~isfield(opt, 'allowDehydrogenaseKnockout')
        opt.allowDehydrogenaseKnockout = true;
    end

    % name global variables
    global maxTime printIntermediateSolutions
    global intermediateSolutionsFile useCobraSolver
    
    % set local variables
    knockType = opt.knockType;
    targetRxn = opt.targetRxn;
    swapNum = opt.swapNum;
    knockoutNum = opt.knockoutNum;
    interventionNum = opt.interventionNum;
    knockableRxns = opt.knockableRxns;
    notKnockableRxns = opt.notKnockableRxns;
    useCobraSolver = opt.useCobraSolver;
    maxW = opt.maxW;
    biomassRxn = opt.biomassRxn;
    dhRxns = opt.dhRxns;
    maxTime = opt.maxTime;
    printIntermediateSolutions = opt.printIntermediateSolutions;
    intermediateSolutionsFile = opt.intermediateSolutionsFile;
    allowDehydrogenaseKnockout = opt.allowDehydrogenaseKnockout;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    printLevel = 10;
    
    if (knockType == 2)
        % perform swaps
        [model, newNames, qsCoupling] = modelSwap(model, dhRxns, true);
    else
        qsCoupling = [];
    end

    chemicalInd = find(ismember(model.rxns, targetRxn));
    if printLevel>=3, display(sprintf('chemical index %d', chemicalInd)); end
    
    %parameters
    findMaxWFlag=0;
    P=1;
    coupledFlag = 1;
    K = knockoutNum;
    L = swapNum;
    X = interventionNum;
    
    disp('preparing model')
    consModel = prepareOptSwapModel(model, chemicalInd, biomassRxn);

    disp('createRobustOneWayFluxesModel')
    [model, yInd, yCoupledInd, qInd, qCoupledInd, sInd, sCoupledInd,...
     notYqsInd, notYqsCoupedInd, m, n, coupled] = ...
        createRobustOneWayFluxesModel(consModel, chemicalInd, coupledFlag, ...
                                      knockableRxns, notKnockableRxns, ...
                                      qsCoupling, useCobraSolver);
  
    
    % sizes
    ySize = length(yInd); yCoupledSize = length(yCoupledInd);
    qSize = length(qInd); qCoupledSize = length(qCoupledInd);
    sSize = length(sInd); sCoupledSize = length(sCoupledInd);
    if (qSize ~= sSize)
        error('Dehydrogenases do not match swap reactions.');
    end 
    if printLevel>=3
        display(sprintf('ySize %d, qSize %d, sSize %d',...
                        ySize, qSize, sSize)); 
        display(sprintf('yCoupledSize %d, qCoupledSize %d, sCoupledSize %d',...
                        yCoupledSize, qCoupledSize, sCoupledSize)); 
    end 
    
    % don't sort combination indexes
    % yqsCoupledInd = [yCoupledInd; qCoupledInd; sCoupledInd];
    yqsInd = [yInd; qInd; sInd];
    
    % combined sizes
    yqsSize = ySize + qSize + sSize; 
    yqsCoupledSize = yCoupledSize + qCoupledSize + sCoupledSize;
    
    % optSwap diverges here
    if (knockType==0 || knockType==1)

        coupledYSize = size(yCoupledInd, 1);
        coupledYInd = yCoupledInd;
        notYInd = notYqsInd; coupledNotYInd = notYqsCoupedInd;
        
    
        %part 2
        %max C'v
        %s.t
        %[A,Ay]*[v;y]<=B
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
        selYInd=zeros(m,1);
        selYInd(yInd)=1;

        Ay1=diag(selYInd);
        Ay1(coupled(:,2), :)=Ay1(coupled(:,1), :);
        Ay1=Ay1*diag(model.ub);
        for j=1:length(coupled)
            if (model.ub(coupled(j,1)) ~=0)
                Ay1(coupled(j,2), coupled(j,1))=Ay1(coupled(j,2), coupled(j,1)).*(model.ub(coupled(j,2))./model.ub(coupled(j,1)));
            end
        end

        Ay2=diag(selYInd);
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

        Ay=[zeros(2*n+2*(vSize-ySize-coupledYSize),ySize);
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
           zeros(2*(ySize+coupledYSize),1);    ];

        C=model.organismObjective;

        %needs to be minimum so:
        %min -C*v
        %s.t
        %-[A,Ay]*[v;y]>=-B
        disp('separateTransposeJoin')
        [A_w, Ay_w ,B_w,C_w, lb_w, ub_w, wSize, wZs] = ...
            separateTransposeJoin(-A, -Ay,-B,-C,ySize, 1,  m, maxW,findMaxWFlag, zSize);
        %max C_w(w  z)'
        %s.t
        %[A_w  Ay_w]*(w z y)  <=  B_w
        awSizeRow=size(A_w,1);

        Ajoined=[
        %C', P*C_w', sparse(1, ySize);    %dual problem objective function = primal problem objective function
            -C', -P*C_w', sparse(1, ySize);
            A, sparse(aSizeRow, wSize+zSize), Ay;     %stochiometric constraints + flux boundry constraints
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



        if (knockType == 0)
            disp('optKnock')

            results = setupAndRunMILP(Cjoined, Ajoined, Bjoined, ...
                                      lbJoined, ubJoined, IntVars_optKnock, ...
                                      model, yInd, [], [], K, [], ...
                                      coupledFlag, coupled);

        elseif (knockType == 1)
            disp('robustKnock')
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

            disp('separateTransposeJoin')
            [A2_w, Ay2_w ,B2_w,C2_w, lb2_w, ub2_w, uSize, uZs]=...
                separateTransposeJoin(A2, Ay2,B2,C2 ,ySize, 1,  vSize+wSize+zSize,...
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

            results = setupAndRunMILP(C3, A3, B3, lb3, ub3, intVars,...
                                      model, yInd, [], [], K, [], ...
                                      coupledFlag, coupled); 
        end
        
    elseif (knockType == 2) 
        disp('optSwap');
        
        %part 2
        %max C'v
        %s.t
        %[A,Ay]*[v;y]<=B
        I = eye(m);
        A=[ 
            model.S;
            -model.S;
            I(notYqsInd,:);
            I(notYqsCoupedInd,:);
            -I(notYqsInd,:);
            -I(notYqsCoupedInd,:);
            I([yInd; yCoupledInd; qInd; qCoupledInd; sInd; sCoupledInd;],:); 
            -I([yInd; yCoupledInd; qInd; qCoupledInd; sInd; sCoupledInd;],:);
          ];

        [aSizeRow, vSize] = size(A);
        selYInd = zeros(m,1); selQInd = zeros(m,1); selSInd = zeros(m,1);
        selYInd(yInd) = 1; selQInd(qInd) = 1; selSInd(sInd) = 1;
        
        % Ay1, Aq1, As1
        Ay1=diag(selYInd);
        Ay1(coupled(:,2), :)=Ay1(coupled(:,1), :);
        Ay1=Ay1*diag(model.ub);     
        
        Aq1=diag(selQInd);
        Aq1(coupled(:,2), :)=Aq1(coupled(:,1), :);
        Aq1=Aq1*diag(model.ub);
        
        As1=diag(selSInd);
        As1(coupled(:,2), :)=As1(coupled(:,1), :);
        As1=As1*diag(model.ub);
        
        for j=1:length(coupled)
            if (model.ub(coupled(j,1)) ~=0)
                Ay1(coupled(j,2), coupled(j,1)) = ...
                           Ay1(coupled(j,2), coupled(j,1))  .*  ...
                           (model.ub(coupled(j,2)) ./ model.ub(coupled(j,1)));
                Aq1(coupled(j,2), coupled(j,1)) = ...
                           Aq1(coupled(j,2), coupled(j,1))  .*  ...
                           (model.ub(coupled(j,2)) ./ model.ub(coupled(j,1)));
                As1(coupled(j,2), coupled(j,1)) = ...
                           As1(coupled(j,2), coupled(j,1))  .*  ...
                           (model.ub(coupled(j,2)) ./ model.ub(coupled(j,1)));
            end
        end
        
        % Ay2, Aq2, As2
        Ay2=diag(selYInd);
        Ay2(coupled(:,2), :)=Ay2(coupled(:,1), :);
        Ay2=Ay2*diag(model.lb);     
        
        Aq2=diag(selQInd);
        Aq2(coupled(:,2), :)=Aq2(coupled(:,1), :);
        Aq2=Aq2*diag(model.lb);
        
        As2=diag(selSInd);
        As2(coupled(:,2), :)=As2(coupled(:,1), :);
        As2=As2*diag(model.lb);
        
        for j=1:length(coupled)
            if (model.lb(coupled(j,1)) ~=0)
                Ay2(coupled(j,2), coupled(j,1)) = ...
                           Ay2(coupled(j,2), coupled(j,1))  .*  ...
                           (model.lb(coupled(j,2)) ./ model.lb(coupled(j,1)));
                Aq2(coupled(j,2), coupled(j,1)) = ...
                           Aq2(coupled(j,2), coupled(j,1))  .*  ...
                           (model.lb(coupled(j,2)) ./ model.lb(coupled(j,1)));
                As2(coupled(j,2), coupled(j,1)) = ...
                           As2(coupled(j,2), coupled(j,1))  .*  ...
                           (model.lb(coupled(j,2)) ./ model.lb(coupled(j,1)));
            end
        end

        z1 = [find(Ay1); find(Aq1); find(As1)];
        z2 = [find(Ay2); find(Aq2); find(As2)];
        zSize = size([z1;z2],1);
        
        if printLevel>=3, disp(sprintf('zSize %d', zSize)); end
        
        yqsCoupledSize = length(yCoupledInd) + length(qCoupledInd) + length(sCoupledInd);
        Ayqs = [
            zeros(2*n+2*(vSize-yqsSize-yqsCoupledSize),yqsSize); 
            -Ay1(yInd,yqsInd);
            -Ay1(yCoupledInd,yqsInd);  
            -Aq1(qInd,yqsInd);
            -Aq1(qCoupledInd,yqsInd);
            -As1(sInd,yqsInd);
            -As1(sCoupledInd,yqsInd);
            Ay2(yInd,yqsInd);
            Ay2(yCoupledInd,yqsInd);     
            Aq2(qInd,yqsInd);
            Aq2(qCoupledInd,yqsInd);
            As2(sInd,yqsInd);
            As2(sCoupledInd,yqsInd);
             ];  %flux boundary constraints
         
        
        %so: [A,Ayqs]x<=B;
        B = [
            zeros(2*n,1);
           model.ub(notYqsInd);
           model.ub(notYqsCoupedInd);
           -model.lb(notYqsInd);
           -model.lb(notYqsCoupedInd);
            zeros(2 * (yqsSize + yqsCoupledSize), 1);    
          ];

        C=model.organismObjective;

        %needs to be minimum so:
        %min -C*v
        %s.t
        %-[A,Ayqs]*[v;y]>=-B
        disp('separateTransposeJoin')
        [A_w, Ayqs_w ,B_w,C_w, lb_w, ub_w, wSize, wZs] = ...
            separateTransposeJoin(-A, -Ayqs, -B, -C, yqsSize, ...
                                  1,  m, maxW, findMaxWFlag, zSize);
        %max C_w(w  z)'
        %s.t
        %[A_w  Ayqs_w]*(w z y)  <=  B_w
        
        awSizeRow = size(A_w, 1);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%
        
        

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

        A2 = - [
            -C', -P*C_w';
            A, sparse(aSizeRow, wSize + zSize);
            sparse(awSizeRow, vSize), A_w
               ];

        Ayqs2 = - [
            zeros(1, yqsSize);
            Ayqs;
            Ayqs_w
               ];

        C2 = [
            model.C_chemical;
            zeros(wSize + zSize,1)
                  ];
        
        B2 = - [
            0;
            B;
            B_w;
               ];

        z3 = find(Ayqs2);
        zSizeOptKnock2=size(z3,1);

        disp('separateTransposeJoin')
        [A2_w, Ayqs2_w ,B2_w,C2_w, lb2_w, ub2_w, uSize, ~]=...
            separateTransposeJoin(A2, Ayqs2, B2, C2, yqsSize, ...
                                  1,  vSize+wSize+zSize, ...
                                  maxW,findMaxWFlag, zSizeOptKnock2); 
        
        %max C2_w*x'
        %s.t
        %A2_w*x+Ayqs2_w*y  <=  B2_w

        %add u1, u2 variables so y will be feasible. add the constraints:
        % su=0 and umin*y<u<umax*y
        [A2_wRow, A2_wCol]=size(A2_w);
        [ARow, ACol] = size(A);

        sCoupledMatrix = zeros(qSize, qSize);
        if ~isempty(qsCoupling)
            qsCoupling_S = qsCoupling(:,2);
        end
        for i = 1:qSize
            thisSIndex = qsCoupling_S(qsCoupling(:,1)==qInd(i));
            sCoupledMatrix(i,sInd==thisSIndex) = 1;
        end
        A3 = [
        %dual constraints
            A2_w, sparse(A2_wRow,ACol), Ayqs2_w;   %  19316       23762
        % swap constraint
            zeros(qSize, uSize + zSizeOptKnock2 + vSize + ySize), eye(qSize), sCoupledMatrix;
        %feasibility constraint
            sparse(ARow,uSize+zSizeOptKnock2), A, Ayqs   % 10134       23762
           ];

        B3=[
        % dual constraints
            B2_w;
        % swap constraint
            ones(qSize, 1);
        % feasibility constraint
            B
           ];
        
        % require swap if property is false
        if ~allowDehydrogenaseKnockout
            display('dehydrogenase knockouts not allowed');
             A3 = [A3;
                  zeros(qSize, uSize + zSizeOptKnock2 + vSize + ySize), ...
                  -eye(qSize), -sCoupledMatrix
                 ];
             B3 = [B3; -ones(qSize, 1)];
        end
        
        % add knockout and swap count constraints
        if K >= 0 % knockoutNum
            A3 = [A3; ...
                  zeros(1, uSize + zSizeOptKnock2 + vSize), ...
                  -ones(1, ySize + qSize + sSize)
                 ];
            B3 = [B3; K - ySize - qSize];
        end
        if L >= 0 % swapNum
            A3 = [A3; ...
                  zeros(1, uSize + zSizeOptKnock2 + vSize + ySize), ...
                  zeros(1, qSize),...
                  ones(1, sSize)];
            B3 = [B3; L];
        end
        if X >= 0 % interventionNum
            A3 = [A3; ...
                  zeros(1, uSize + zSizeOptKnock2 + vSize), ...
                  -ones(1, ySize + qSize), ...
                  zeros(1, sSize)];
            B3 = [B3; X - ySize - qSize];
        end
        
        C3=[
            C2_w;
            zeros(ACol, 1);
            zeros(yqsSize, 1)
           ];

        tmpL = lb2_w(1:uSize + zSizeOptKnock2);
        tmpH=ub2_w(1:uSize + zSizeOptKnock2);
        ysUpperBound = ones(yqsSize, 1);
        lb3 = [
            tmpL; 
            model.lb; 
            zeros(yqsSize,1)
              ];
        ub3 = [
            tmpH; 
            model.ub; 
            ysUpperBound
              ];
        intVars = (A2_wCol + ACol + 1):(A2_wCol + ACol + yqsSize);
        
        results = setupAndRunMILP(C3, A3, B3, lb3, ub3, intVars, ...
                                  model, yInd, qInd, sInd, K, L, ...
                                  coupledFlag, coupled);
    end
end

function results = setupAndRunMILP(C, A, B, lb, ub, intVars,...
                                   model, yInd, qInd, sInd, K, L, ...
                                   coupledFlag, coupled);

    global maxTime
    global printIntermediateSolutions intermediateSolutionsFile
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
    
    ySize = length(yInd); qSize = length(qInd); sSize = length(sInd);

    %solving mip
    disp('set up MILP')
    global useCobraSolver
    if (useCobraSolver)
        % maximize
        MILPproblem.c = C;
        MILPproblem.osense = -1; % max
        MILPproblem.A = A;
        MILPproblem.b_L = [];
        MILPproblem.b_U = B;
        MILPproblem.b = B;
        MILPproblem.lb = lb; MILPproblem.ub = ub;
        MILPproblem.x0 = [];
        MILPproblem.vartype = char(ones(1,length(C)).*double('C'));
        MILPproblem.vartype(intVars) = 'I';
        MILPproblem.csense = char(ones(1,length(B)).*double('L'));
        
        if printIntermediateSolutions
            disp('setting up callback')
            warning('not finished')
        end
        [MILPproblem, solverParams] = setParams(MILPproblem, true, maxTime); 
        disp('Run COBRA MILP')
        Result_cobra = solveCobraMILP(MILPproblem,solverParams);

        results.raw = Result_cobra;
        results.y = Result_cobra.full(intVars(1:ySize));
        results.q = Result_cobra.full(intVars(ySize+1:ySize+qSize));
        results.s = Result_cobra.full(intVars(ySize+qSize+1:ySize+qSize+sSize));
        results.exitFlag = Result_cobra.stat;
        results.inform = Result_cobra.origStat;
        results.chemical = -Result_cobra.obj;
        results.f_k = Result_cobra.obj;
        results.solver = Result_cobra.solver;
        save results
    else
        % minimize
        Prob_OptKnock2=mipAssign(-C, A, [], B, lb, ub, [], 'part 3 MILP', ...
                                 setupFile, nProblem, ...
                                 intVars, VarWeight, KNAPSACK, fIP, xIP, ...
                                 f_Low, x_min, x_max, f_opt, x_opt);
        disp('setParams')
        Prob_OptKnock2 = setParams(Prob_OptKnock2, false, maxTime);

        if printIntermediateSolutions
            disp('setting up callback')
            global optSwapCallbackOptions
            optSwapCallbackOptions.outputFile = intermediateSolutionsFile; 
            optSwapCallbackOptions.model = model;
            optSwapCallbackOptions.yInd = yInd;
            optSwapCallbackOptions.qInd = qInd;
            optSwapCallbackOptions.intVars = intVars;
            Prob_OptKnock2.MIP.callback(14) = 1;
            Prob_OptKnock2.MIP.callback(7) = 1;
        end
        
        disp('tomRun') 
        warning('hack to show script hierarchy')
        Result_tomRun = tomRun('cplex', Prob_OptKnock2, 10);

        save('raw results','Result_tomRun');

        results.raw = Result_tomRun; 
        results.y = Result_tomRun.x_k(intVars(1:ySize));
        results.q = Result_tomRun.x_k(intVars(ySize+1:ySize+qSize));
        results.s = Result_tomRun.x_k(intVars(ySize+qSize+1:ySize+qSize+sSize));
        results.exitFlag = Result_tomRun.ExitFlag;
        results.inform = Result_tomRun.Inform;
        results.chemical = -Result_tomRun.f_k;
        results.f_k = Result_tomRun.f_k;
        results.solver = 'tomlab_cplex';
    end
    
    results.C = C;
    results.A = A;
    results.B = B;
    results.lb = lb;
    results.ub = ub;
    results.K = K;
    results.L = L; 
    results.intVars = intVars;
    results.yInd = yInd;
    results.qInd = qInd;
    results.sInd = sInd;
    results.organismObjectiveInd = model.organismObjectiveInd;
    results.chemicalInd = model.chemicalInd;
    results.knockoutRxns = model.rxns(yInd(results.y==0));
    results.knockoutDhs = model.rxns(qInd(results.q==0 & results.s==0));
    results.knockoutRxns = [results.knockoutRxns; results.knockoutDhs];
    results.swapRxns = model.rxns(qInd(results.s==1));
    % save('sorted results', 'results');
end



function consModel=prepareOptSwapModel(model, chemicalInd, biomassRxn)

%lets start!

    consModel=model;
    [metNum, rxnNum] = size(consModel.S);
    consModel.row_lb=zeros(metNum,1);
    consModel.row_ub=zeros(metNum,1);
    
    % setup chemical objective
    consModel.C_chemical=zeros(rxnNum,1);
    consModel.C_chemical(chemicalInd)=1;
    
    % add chemical index to model
    consModel.chemicalInd = chemicalInd;

    % setup biomass objective
    biomassInd = find(ismember(model.rxns, biomassRxn));
    if isempty(biomassInd)
        error('biomass reaction not found in model');
        return;
    end
    consModel.organismObjectiveInd = biomassInd;
    consModel.organismObjective = zeros(rxnNum,1);
    consModel.organismObjective(biomassInd) = 1;
    
    %remove small  reactions
    sel1 = (consModel.S>-(10^-3) & consModel.S<0);
    sel2 = (consModel.S<(10^-3) & consModel.S>0);
    consModel.S(sel1)=0;
    consModel.S(sel2)=0;

end