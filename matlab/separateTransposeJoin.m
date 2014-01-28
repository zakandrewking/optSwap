function [A_w,Ay_w ,B_w,C_w, lb_w, ub_w, wSize, wZs]=separateTransposeJoin(A, Ay, B,C,ySize, yMax, aSizeCol, highNum, findMaxFlag, zSize)

%Return maximum function under the following assumptions:
%1. y variables are the last variables at A (the last columns).
%2. y is not in the objective function.
%3. each original value v has a relation with one y variable at the most
%4. it is a minimum problem with constraint Ax>=B

maxW=0;
%separation
%min c1*v
%s.t
%A*v >= B-Ay*y

%transpose
%max ( B-Ay*y)'w=B'w-(Ay*y)'w
%s.t
%A'*w<=c
wSize=size(A',2);

%force linearity by: z=y*w
%max ( B-Ay*y)'w=B'w-cy'z
%s.t
%A'*w<=c
%D(w,z)<=H

wMax=highNum*ones(wSize,1);

HConst=zeros(3*zSize,1);
Hy=zeros(3*zSize,ySize);
D=sparse(3*zSize,wSize+zSize);
Cz=zeros(zSize,1);
zInd=1;
wZs=zeros(zSize,1);

d1=[0;1;-1];
d2=[1;-1;1];
hc=[0;1;0];
hy=[-1;1;0];

for k=1:wSize
    [maxVal,indMax]=max(Ay(k,:));
    [minVal,indMin]=min(Ay(k,:));
    if (maxVal~=0)
        val=maxVal;
        ind=indMax;
    end
    if (minVal~=0)  %for minus values
        val=minVal;
        ind=indMin;
    end
    if (maxVal~=0  || minVal~=0)
        if (findMaxFlag==1)
            c_k=zeros(wSize,1); %find maximal value
            c_k(k)=1;
            Prob=lpAssign(-c_k, A', [], C, zeros(wSize,1), []);
            Result=tomRun('cplex', Prob, 0);
            exitFlag=Result.ExitFlag;
            res = -Result.f_k;
            if (exitFlag==0 && res > 1e-3) %found a solution
                wMax(k)=res;
                if (res > maxW)
                    maxW=res;
                end
            end
        end
        constInd=3*zInd-2:3*zInd;
        Cz(zInd,1)=val;
        D(constInd,k)=d1;
        D(constInd,zInd+wSize)=d2;
        HConst(constInd)=wMax(k,1)*hc;
        Hy(constInd,ind)=wMax(k,1)*hy;
        wZs(zInd)=k;
        zInd=zInd+1;
    end
end

%new objective function
%the new variables are (w,z,y)
%max C_w(w  z)'
%s.t
%[A_w  Ay_w]*(w z y)  <=  B_w

C_w=[B;
    -Cz
    ];
A_w=[
    A', zeros(aSizeCol, zSize);
    D;
    ];
Ay_w=[
    zeros(aSizeCol, ySize);
    Hy;
    ];
B_w=[
    C;
    HConst;
    ];
lb_w=zeros(wSize+zSize+ySize, 1);
ub_w=[highNum*ones(wSize+zSize, 1); yMax*ones(ySize,1)];


