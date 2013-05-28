% log likelihood ratio function, with non-linear equality and inequality constraits
% April 29, 2013

function [LLR,a,varargout]=LogLikelihoodRatioNonc2(x,inParams)
%function LLR=LogLikelihoodRatio(x,inParams)
xmaxmin = inParams.xmaxmin;

%global Np alphaP deltaP kp N timingResiduals sd yr
%global variables defined in the simulator: alphaP(i),deltaP(i),yr,s
%disp(['number of elements in x: ', num2str(numel(x))]);
% 0<=x<=1, convert to the physical range of each variable
% alpha=x(1)*(xmaxmin(1,1)-xmaxmin(1,2))+xmaxmin(1,2); % *2*pi;  % [0, 2*pi]
% delta=x(2)*(xmaxmin(2,1)-xmaxmin(2,2))+xmaxmin(2,2);    %*pi-pi/2.0;  % [-pi/2, pi/2]
% omega=x(3)*(xmaxmin(3,1)-xmaxmin(3,2))+xmaxmin(3,2);   %*10+10;  % [10, 20]  ??
% phi0=x(4)*(xmaxmin(4,1)-xmaxmin(4,2))+xmaxmin(4,2);  %*2*pi;  % [0, 2*pi]

% transform the standard coordinates [0,1] back to physical coordinates
minAlpha = xmaxmin(1,2);
rangeAlpha = xmaxmin(1,1)-minAlpha;
minDelta = xmaxmin(2,2);
rangeDelta = xmaxmin(2,1)-minDelta;
alpha=x(1)*rangeAlpha+minAlpha;  % [0, 2*pi]
delta=x(2)*rangeDelta+minDelta;  % [-pi/2, pi/2]

omega=x(3)*(xmaxmin(3,1)-xmaxmin(3,2))+xmaxmin(3,2);  % [10, 20]  ??
phi0=x(4)*(xmaxmin(4,1)-xmaxmin(4,2))+xmaxmin(4,2);

phiI = zeros(1,8);  % arbitrary phase
phiI(1)=x(5)*(xmaxmin(5,1)-xmaxmin(5,2))+xmaxmin(5,2);
phiI(2)=x(6)*(xmaxmin(6,1)-xmaxmin(6,2))+xmaxmin(6,2);
phiI(3)=x(7)*(xmaxmin(7,1)-xmaxmin(7,2))+xmaxmin(7,2);
phiI(4)=x(8)*(xmaxmin(8,1)-xmaxmin(8,2))+xmaxmin(8,2);
phiI(5)=x(9)*(xmaxmin(9,1)-xmaxmin(9,2))+xmaxmin(9,2);
phiI(6)=x(10)*(xmaxmin(10,1)-xmaxmin(10,2))+xmaxmin(10,2);
phiI(7)=x(11)*(xmaxmin(11,1)-xmaxmin(11,2))+xmaxmin(11,2);
phiI(8)=x(12)*(xmaxmin(12,1)-xmaxmin(12,2))+xmaxmin(12,2);

% phiI = zeros(1,8);  % arbitrary phase
% phiI(1)=x(5)*2*pi;
% phiI(2)=x(6)*2*pi;
% phiI(3)=x(7)*2*pi;
% phiI(4)=x(8)*2*pi;
% phiI(5)=x(9)*2*pi;
% phiI(6)=x(10)*2*pi;
% phiI(7)=x(11)*2*pi;
% phiI(8)=x(12)*2*pi;

%if nargin == 3
%     if alpha > xmaxmin(1,1) || alpha < xmaxmin(1,2)
%        LLR=Inf;
%        if nargout > 1
%            varargout{1}=[alpha,delta,omega,phi0,phiI];
%        end
%        return
%     end
%     
%     if delta > xmaxmin(2,1) || delta < xmaxmin(2,2)
%        LLR=Inf;
%        if nargout > 1
%            varargout{1}=[alpha,delta,omega,phi0,phiI];
%        end
%        return
%     end
%end


if nargout > 1
    varargout{1}=[alpha,delta,omega,phi0,phiI];
end

% for structure inParams
Np = inParams.Np;
N = inParams.N;
s = inParams.s;
sd = inParams.sd;
alphaP = inParams.alphaP;
deltaP = inParams.deltaP;
%kp = inParams.kp;
yr = inParams.yr;

N1=zeros(8,1);  % 1-D vector
M2=zeros(8,8);  % 2-D matrix
Abase=zeros(Np,8,N);  % Np*8*N cube, base functions for each pulsar


for i=1:1:Np
    Abase(i,:,:)=AbaseFunc(alpha,delta,omega,phi0,phiI(i),alphaP(i),deltaP(i),yr);
end

tmp=zeros(1,N);
for mu=1:1:8

    % calculate 8*1 vector N
    for i=1:1:Np
        tmp(1,:)=Abase(i,mu,:);  % convert (1,1,N) matrix to (1,N)
        %size(tmp)
        N1(mu)=N1(mu)+InnProduct(s(i,:),tmp,sd);
        %N1(mu)=N1(mu)+InnProduct(tmp1,tmp2,sd);
    end
    
    % calculate 8*8 matrix M
    for nu=1:1:8
        
        for i=1:1:Np
            M2(mu,nu)=M2(mu,nu)+InnProduct(Abase(i,mu,:),Abase(i,nu,:),sd);
        end
        
    end

end


B=diag([1,-1,0,0,1,-1,0,0]);
% C is the constraint matrix
C=[1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0; ...
   0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0; ...
   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0; ...
   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0];

% inverse square matrix, warning message is printed if X is badly scaled or nearly singular.
a1=M2\N1;  % tmp variable, \ backslash or matrix left division, more efficient
a2=M2\C';
% a: Lagrangian multiplier solution with only 4 linear equality constraints, 
% ref eq.5 of my note on 4/23/13
a=a1-a2*((C*a2)\(C*a1));
abest=a;
fbest = -0.5 * ( dot(N1,a1) - (C*a1)'*((C*a2)\(C*a1)) );
%M2inv=zeros(8,8);
%M2inv=inv(M2);

% LLR with 4 linear constraints, and additional nonlinear constraints are 
% subject if certain condition is violated. Here accually LLR = - LLR

% setup for local optimization 'fmincon'
%a0 = a;  % make a starting guess at the solution
%beq=[0.0; 0.0; 0.0; 0.0];  % right hand side of 4 linear constraints
%options = optimset('Algorithm','active-set','Display','off');
options = optimset('Algorithm','interior-point','Display','off'); %, ...
    %'Hessian','user-supplied','HessFcn',@hessianfcn);
%func = @(x) exp(x(1))*(4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 2*x(2) + 1);
%func = @(a) -(a'*N1-0.5*a'*M2*a);  % min log likelihood ratio
Ampfunc = @ (a) 0.5*(sqrt(a(1)^2+a(5)^2)+sqrt(a(1)^2+a(5)^2-a(2)^2-a(6)^2));
a0=zeros(4,1); % initial value for nonlinear constraint quadratic problem
penalty=0.9;  % penalty factor when every trial is failed


if a'*B*a > 0
    Amp=Ampfunc(a);
    
    if a(2)^2-4*Amp^2*a(5)^2/(a(1)^2+a(5)^2) <= 0
        % CASE 1: linear constraits and both nonlinear constraints are satisfied
        
        %LLR=-0.5*dot(N1,a1);  % LLR without constraints
        %LLR = -0.5 * ( N1'*M2inv*N1 - (C*M2inv*N1)'*(C*M2inv*C')*(C*M2inv*N1) );
        %LLR = -0.5 * ( dot(N1,a1) - (C*a1)'*inv(C*M2inv*C')*(C*a1) );
        LLR = fbest; % -0.5 * ( dot(N1,a1) - (C*a1)'*((C*a2)\(C*a1)) );  % maybe most efficient and accurate
        %LLR = -0.5 * ( dot(N1,a1) - (C*a1)'*inv(C*a2)*(C*a1) );
        %LLR = -0.5 * ( dot(N1,a1) - (C*a1)'*(C*a2)\(C*a1) );
    else
        % CASE 2: second nonlinear constraint is violated
        %a0=a; %+a.*(0.1*rand(8,1)-0.05);  % make a starting guess at the solution
        a0(1)=a(1);
        a0(2)=a(2);
        a0(3)=a(5);
        a0(4)=a(6);
        
        %[a,fval] = fmincon(func,a0,[],[],C,beq,[],[],@nonlconst2,options);
        [y,fval] = fmincon(@(x) func4(x,N1,M2),a0,[],[],[],[],[],[],@nonlconst2,options);
        %[a,fval] = fmincon(func,a0,[],[],C,beq,[],[],@nonlconst1,options);
        a(1)=y(1);
        a(2)=y(2);
        a(3)=-y(1);
        a(4)=y(2);
        a(5)=y(3);
        a(6)=y(4);
        a(7)=-y(3);
        a(8)=y(4);
        
        if a'*B*a >= 0 && a(2)^2-4*Ampfunc(a)^2*a(5)^2/(a(1)^2+a(5)^2) <= 0
            LLR=fval;
        else
            %LLR=Inf;  % may be not able to get a decent solution
            %LLR=0.0;
            a=abest;
            LLR=fbest*penalty; %fval;
            %disp('In LoglikelihoodRatioNonc.m: Inf occurs at 1');
        end
        
    end
    
else
    
    %disp('Warning: In LogLikelihoodRatio.m: a(1)^2+a(5)^2-a(2)^2-a(6)^2 < 0');
    % CASE 3:
    %a0=a+a.*(0.1*rand(8,1)-0.05);
    a0(1)=a(1);
    a0(2)=a(2);
    a0(3)=a(5);
    a0(4)=a(6);

    %[a,fval] = fmincon(func,a0,[],[],C,beq,[],[],@nonlconst1,options);
    [y,fval] = fmincon(@(x) func4(x,N1,M2),a0,[],[],[],[],[],[],@nonlconst1,options);
    a(1)=y(1);
    a(2)=y(2);
    a(3)=-y(1);
    a(4)=y(2);
    a(5)=y(3);
    a(6)=y(4);
    a(7)=-y(3);
    a(8)=y(4);
    
    % fmincon may be not accurate, ...
    
    if a'*B*a >= 0 && a(2)^2-4*Ampfunc(a)^2*a(5)^2/(a(1)^2+a(5)^2) <= 0 % >
        LLR=fval;
    else
        %a0=a+a.*(0.1*rand(8,1)-0.05);
        a0(1)=a(1);
        a0(2)=a(2);
        a0(3)=a(5);
        a0(4)=a(6);
        
        %[a,fval] = fmincon(func,a0,[],[],C,beq,[],[],@nonlconst2,options);
        %[a,fval] = fmincon(func,a0,[],[],C,beq,[],[],@nonlconst1,options);
        [y,fval] = fmincon(@(x) func4(x,N1,M2),a0,[],[],[],[],[],[],@nonlconst2,options);
        a(1)=y(1);
        a(2)=y(2);
        a(3)=-y(1);
        a(4)=y(2);
        a(5)=y(3);
        a(6)=y(4);
        a(7)=-y(3);
        a(8)=y(4);
        
        if a'*B*a >= 0 && a(2)^2-4*Ampfunc(a)^2*a(5)^2/(a(1)^2+a(5)^2) <= 0
            LLR=fval;
        else
            %LLR=Inf;  % may be not able to get a decent solution
            %LLR=0.0;
            a=abest;
            LLR=fbest*penalty;  %fval;
            %disp('In LoglikelihoodRatioNonc.m: Inf occurs at 2');

        end
        
    end
    
end

% END of function