function results = optSwapD(model, opt)
% Adapted from RobustKnock.m
% By Zachary King 8/13/2012
%
% INPUTS
% model
% dhRxns
%
% OPTIONAL
% opt - struct with an of the following options:
%   knockType - 2 run OptSwap; 1 run robustKnock; 0 run OptKnock (broken)
%   targetRxn - chemical to produce
%   knockoutNum - number of knockouts
%   knockableRxns - reactions that should definitely be considered for knockout
%   notKnockableRxns - reactions that should be excluded from knockout
%   useCobraSolver - (unfinished) 1 use any cobra solver; 0 use tomlab cplex solver
%   maxW - maximal value of dual variables (higher number will be more
%          accurate but takes more calculation time)
%   biomassRxn 
%   % FOR OPTSWAP:
%   swapNum - number of swaps
%   dhRxns - dehydrogenase reaction list
%
% OUTPUTS
% results


    if nargin < 1, error('Not enough arguments'); end

    % set default parameters
    if ~exist('opt','var')
        opt = struct();
    end
    if ~isfield(opt,'targetRxn'), opt.targetRxn = 'EX_for(e)'; end
    if ~isfield(opt,'knockoutNum'), opt.knockoutNum = 3; end
    if ~isfield(opt,'knockType'), opt.knockType = 0; end
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
    if ~isfield(opt,'swapNum'), opt.swapNum = 0; end
    if ~isfield(opt,'dhRxns'), opt.dhRxns = {}; end
    
    printLevel = 10;
    
    if (opt.knockType == 2)
        % perform swaps
        [model, newNames, qsCoupling] = modelSwap(model, opt.dhRxns, true);
    else
        qsCoupling = [];
    end

    chemicalInd = find(ismember(model.rxns, opt.targetRxn));
    if printLevel>=3, display(sprintf('chemical index %d', chemicalInd)); end
    
    %parameters
    findMaxWFlag=0;
    P=1;
    coupledFlag = 1;
    K = opt.knockoutNum;
    L = opt.swapNum;
    
    disp('preparing model')
    consModel = prepareOptSwapModel(model, chemicalInd, opt.biomassRxn);

    disp('createRobustOneWayFluxesModel')
    [model, yInd, yCoupledInd, qInd, qCoupledInd, sInd, sCoupledInd,...
     notYqsInd, notYqsCoupedInd, m, n, coupled] = ...
                createRobustOneWayFluxesModel(consModel, chemicalInd, coupledFlag, ...
                                      opt.knockableRxns, opt.notKnockableRxns, qsCoupling,...
                                              opt.useCobraSolver);
  
    
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
    if (opt.knockType==0 || opt.knockType==1)

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
            separateTransposeJoin(-A, -Ay,-B,-C,ySize, 1,  m, opt.maxW,findMaxWFlag, zSize);
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



        if (opt.knockType == 0)
            disp('optKnock')


            results = setupAndRunMILP(opt.useCobraSolver, ...
                                      Cjoined, Ajoined, Bjoined, ...
                                      lbJoined, ubJoined, IntVars_optKnock, ...
                                      model, yInd, [], [], K, [], ...
                                      coupledFlag, coupled);
            % disp('mipAssign')
            % Prob_optKnock=mipAssign(-Cjoined, Ajoined, [], Bjoined, lbJoined, ubJoined, [], 'part 2.3 MILP', ...
            %                         setupFile, nProblem, ...
            %                         IntVars_optKnock, VarWeight, KNAPSACK, fIP, xIP, ...
            %                         f_Low, x_min, x_max, f_opt, x_opt);
            % disp('setParams')
            % Prob_optKnock=setParams(Prob_optKnock);
            % disp('tomRun')
            % Result_optKnock=tomRun('cplex', Prob_optKnock, 1);

            % YsOptKnockNotRound=Result_optKnock.x_k(IntVars_optKnock);
            % YsOptKnock=round(YsOptKnockNotRound);
            % knockedYind=find(1.-YsOptKnock);
            % knockedVoptKnock=yInd(knockedYind);

            % results.Result_optKnock.exitFlag=Result_optKnock.ExitFlag;
            % results.Result_optKnock.organismObjectiveInd=model.organismObjectiveInd;
            % results.Result_optKnock.chemicalInd=model.chemicalInd;
            % results.Result_optKnock.objective=Result_optKnock.x_k(model.organismObjectiveInd);
            % results.Result_optKnock.chemicel=Result_optKnock.x_k(model.chemicalInd);
            % results.Result_optKnock.x_k=Result_optKnock.x_k;
            % results.Result_optKnock.ySize=ySize;
            % results.Result_optKnock.ySum=sum(YsOptKnockNotRound);
            % results.Result_optKnock.yBoundry=ySize-K;
            % results.Result_optKnock.knockedYind=knockedYind;
            % results.Result_optKnock.knockedYvals=YsOptKnockNotRound(knockedYind);
            % results.Result_optKnock.knockedV=knockedVoptKnock;
            % results.Result_optKnock.knockedVvalues=Result_optKnock.x_k(yInd(results.Result_optKnock.knockedYind));
            % results.Result_optKnock.y=YsOptKnockNotRound;

            % if(Result_optKnock.ExitFlag ==0)        %a sucessful run
            %     disp('optknock successful run')
            %     results.Result_optKnock.knockCheck = ...
            %         knockCheck(model, knockedVoptKnock, coupledFlag, coupled);
            % end




        elseif (opt.knockType == 1)
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
                                      opt.maxW,findMaxWFlag, zSizeOptKnock2);

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

            results = setupAndRunMILP(opt.useCobraSolver, ...
                                      C3, A3, B3, lb3, ub3, intVars,...
                                      model, yInd, [], [], K, [], ...
                                      coupledFlag, coupled);

            % %solving mip
            % disp('mipAssign')
            % Prob_OptKnock2=mipAssign(-C3, A3, [], B3, lb3, ub3, [], 'part 3 MILP', ...
            %                          setupFile, nProblem, ...
            %                          intVars, VarWeight, KNAPSACK, fIP, xIP, ...
            %                          f_Low, x_min, x_max, f_opt, x_opt);
            % disp('setParams')
            % Prob_OptKnock2=setParams(Prob_OptKnock2);
            % disp('tomRun')

            % Result_OptKnock2=tomRun('cplex', Prob_OptKnock2, 2);

            % YsNotRound=Result_OptKnock2.x_k(intVars);
            % Ys=round(YsNotRound);
            % knockedYindNotRound=find(1.-YsNotRound);
            % knockedYindRound=find(1.-Ys);
            % knockedVindRobustKnock=yInd(knockedYindRound);

            % results.Result_OptKnock2.exitFlag=Result_OptKnock2.ExitFlag;
            % results.Result_OptKnock2.organismObjectiveInd=model.organismObjectiveInd;
            % results.Result_OptKnock2.chemicalInd=model.chemicalInd;
            % results.Result_OptKnock2.chemical=-Result_OptKnock2.f_k;
            % results.Result_OptKnock2.exitFlag=Result_OptKnock2.ExitFlag;
            % results.Result_OptKnock2.ySize=ySize;
            % results.Result_OptKnock2.ySum=sum(Result_OptKnock2.x_k(intVars));
            % results.Result_OptKnock2.yLowerBoundry=ySize-K;
            % results.Result_OptKnock2.knockedYindRound=knockedYindRound;
            % results.Result_OptKnock2.knockedYvalsRound=YsNotRound(knockedYindRound);
            % results.Result_OptKnock2.y=Ys;
            % results.Result_OptKnock2.knockedVind=knockedVindRobustKnock;

            % results.Result_OptKnock2.knockCheck = ...
            %     knockCheck(model, knockedVindRobustKnock, coupledFlag, coupled);

        end
        
    elseif (opt.knockType == 2) 
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
        % TODO: where do we use z1, z2?
        if printLevel>=3, disp(sprintf('zSize %d', zSize)); end
        
        yqsCoupledSize = length(yCoupledInd) + length(qCoupledInd) + length(sCoupledInd);
        Ayqs = [
            zeros(2*n+2*(vSize-yqsSize-yqsCoupledSize),yqsSize); 
            -Ay1(yInd,yqsInd);
            -Ay1(yCoupledInd,yqsInd);      % good!! 
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
         
%         yqsIndSort = sort(yqsInd);
%         AyqsSort = [
%             zeros(2*n+2*(vSize-yqsSize-yqsCoupledSize),yqsSize);
%             -Ay1(yInd,yqsIndSort);
%             -Ay1(yCoupledInd,yqsIndSort);      % good!! 
%             -Aq1(qInd,yqsIndSort);
%             -Aq1(qCoupledInd,yqsIndSort);
%             -As1(sInd,yqsIndSort);
%             -As1(sCoupledInd,yqsIndSort);
%             Ay2(yInd,yqsIndSort);
%             Ay2(yCoupledInd,yqsIndSort);     
%             Aq2(qInd,yqsIndSort);
%             Aq2(qCoupledInd,yqsIndSort);
%             As2(sInd,yqsIndSort);
%             As2(sCoupledInd,yqsIndSort);
%              ];  %flux boundry constraints
%          
%         TODO: check this by writing out a few constraint rows
        
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
                                  1,  m, opt.maxW, findMaxWFlag, zSize);
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
                                  opt.maxW,findMaxWFlag, zSizeOptKnock2);

        
        
        
        %max C2_w*x'
        %s.t
        %A2_w*x+Ayqs2_w*y  <=  B2_w

        %add u1, u2 variables so y will be feasible. add the constraints:
        % su=0 and umin*y<u<umax*y
        [A2_wRow, A2_wCol]=size(A2_w);
        [ARow, ACol] = size(A);

        
        A3 = [
        %dual constraints
            A2_w, sparse(A2_wRow,ACol), Ayqs2_w;   %  19316       23762
        %y sum constraints
            zeros(1, uSize + zSizeOptKnock2 + vSize),  -ones(1,ySize), ...
                       zeros(1, qSize + sSize); % 1       23762
        % q sum constraints
            zeros(1, uSize + zSizeOptKnock2 + vSize + ySize), -ones(1, qSize),...
                       zeros(1, sSize); %  1       23762
        % swap constraints
            zeros(qSize, uSize + zSizeOptKnock2 + vSize + ySize), eye(qSize), eye(sSize);
            zeros(qSize, uSize + zSizeOptKnock2 + vSize + ySize), -eye(qSize), -eye(sSize); %  27       23762
        %feasibility constraint
            sparse(ARow,uSize+zSizeOptKnock2), A, Ayqs   % 10134       23762
           ];

        B3=[
            B2_w;
            K - ySize;
            L - qSize; 
            ones(qSize, 1);
            -ones(qSize, 1);
            B
           ];

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
       
        results = setupAndRunMILP(opt.useCobraSolver, ...
                                  C3, A3, B3, lb3, ub3, intVars, ...
                                  model, yInd, qInd, sInd, K, L, ...
                                  coupledFlag, coupled);

    end

end


% fileName2=sprintf('res_%s.mat', datestr(clock, 'yyyymmddTHHMMSS'));
% fl = sprintf('%s//%s', pth, fileName2);
%save (fl);



function results = setupAndRunMILP(useCobraSolver,...
                                   C, A, B, lb, ub, intVars,...
                                   model, yInd, qInd, sInd, K, L, ...
                                   coupledFlag, coupled);

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


    %solving mip
    disp('mipAssign')
    if (useCobraSolver)
        MILPproblem.c = C;
        MILPproblem.osense = -1;
        MILPproblem.A = A;
        MILPproblem.b_L = [];
        MILPproblem.b_U = B;
        MILPproblem.b = B;
        MILPproblem.lb = lb; MILPproblem.ub = ub;
        MILPproblem.x0 = [];
        MILPproblem.vartype = char(ones(1,length(C)).*double('C'));
        MILPproblem.vartype(intVars) = 'I'; % assume no binary 'B'?
        MILPproblem.csense = [];

        [MILPproblem,solverParams] = setParams(MILPproblem, true);

        results = solveCobraMILP(MILPproblem,solverParams);

        % TODO: finish setting up this result analysis:
        YsNotRound=Result_tomRun.x_k(intVars);
        Ys=round(YsNotRound);
        knockedYindNotRound=find(1.-YsNotRound);
        knockedYindRound=find(1.-Ys);
        knockedVindRobustKnock=yInd(knockedYindRound);

        results.exitFlag = Result_tomRun.ExitFlag;
        results.organismObjectiveInd = model.organismObjectiveInd;
        results.chemicalInd = model.chemicalInd;
        results.chemical = -Result_tomRun.f_k;
        results.exitFlag = Result_tomRun.ExitFlag;
        results.ySize = ySize;
        results.ySum = sum(Result_tomRun.x_k(intVars));
        results.yLowerBoundry = ySize-K;
        results.knockedYindRound = knockedYindRound;
        results.knockedYvalsRound = YsNotRound(knockedYindRound);
        results.y = Ys;
        results.knockedVind = knockedVindRobustKnock;

        results.solver = 'cobraSolver';

        % tomlabProblem = mipAssign(osense*c,A,b_L,b_U,lb,ub,x0,'CobraMILP',[],[],intVars);
        % varWeight = []
        % KNAPSACK
        % fIP = []
        % xIP = []
        % f_Low = -1E7
        % x_min = []
        % x_max = []
        % f_opt = -141278
        % x_opt = []

    else

        Prob_OptKnock2=mipAssign(-C, A, [], B, lb, ub, [], 'part 3 MILP', ...
                                 setupFile, nProblem, ...
                                 intVars, VarWeight, KNAPSACK, fIP, xIP, ...
                                 f_Low, x_min, x_max, f_opt, x_opt);

        disp('setParams')
        Prob_OptKnock2 = setParams(Prob_OptKnock2);
        disp('tomRun')

        Result_tomRun = tomRun('cplex', Prob_OptKnock2, 2);

        % save('raw results','Result_tomRun');

        results.raw = Result_tomRun;
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
        
        ySize = length(yInd); qSize = length(qInd); sSize = length(sInd);
        results.y = Result_tomRun.x_k(intVars(1:ySize));
        results.q = Result_tomRun.x_k(intVars(ySize+1:ySize+qSize));
        results.s = Result_tomRun.x_k(intVars(ySize+qSize+1:ySize+qSize+sSize));
        
% TODO: round integer results

%         yqsNotRound=Result_tomRun.x_k(intVars);
%         yqs=round(yqsNotRound);
%         knockedYqsIndNotRound=find(1.-yqsNotRound);
%         knockedYqsIndRound=find(1.-yqs);
%         knockedYqsIndRobustKnock=yqsInd(knockedYqsIndRound);

        results.exitFlag = Result_tomRun.ExitFlag;
        results.inform = Result_tomRun.Inform;
        results.organismObjectiveInd = model.organismObjectiveInd;
        results.chemicalInd = model.chemicalInd;
        results.chemical = -Result_tomRun.f_k;
        results.knockoutRxns = model.rxns(yInd(results.y==0));
        results.knockoutDhs = model.rxns(qInd(results.q==0));
        results.knockoutSwaps = model.rxns(sInd(results.s==0));
%         results.yLowerBoundry = ySize-K;
%         results.knockedYqsIndRound = knockedYqsIndRound;
%         results.knockedYvalsRound = yqsNotRound(knockedYqsIndRound);
%         results.yqs = yqs;
%         results.knockedVind = knockedYqsIndRobustKnock;

%         results.knockCheck = ...
%             knockCheck(model, results.knockedVind, coupledFlag, coupled);

        results.solver = 'tomlab_cplex';
        
        % save('sorted results', 'results');
    end


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