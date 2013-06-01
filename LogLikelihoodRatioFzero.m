% log likelihood ratio function
% April 23, 2013, add constraints for maximizing the likelihood ratio

function [LLR,varargout]=LogLikelihoodRatioFzero(x,inParams)
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
% 
% phiI = zeros(1,8);  % arbitrary phase [0, 2*pi]
% phiI(1)=x(5)*(xmaxmin(5,1)-xmaxmin(5,2))+xmaxmin(4,2);     %2*pi;
% phiI(2)=x(6)*(xmaxmin(6,1)-xmaxmin(6,2))+xmaxmin(4,2);    %2*pi;
% phiI(3)=x(7)*(xmaxmin(7,1)-xmaxmin(7,2))+xmaxmin(4,2);    %2*pi;
% phiI(4)=x(8)*(xmaxmin(8,1)-xmaxmin(8,2))+xmaxmin(4,2);    %2*pi;
% phiI(5)=x(9)*(xmaxmin(9,1)-xmaxmin(9,2))+xmaxmin(4,2);    %2*pi;
% phiI(6)=x(10)*(xmaxmin(10,1)-xmaxmin(10,2))+xmaxmin(4,2);   %2*pi;
% phiI(7)=x(11)*(xmaxmin(11,1)-xmaxmin(11,2))+xmaxmin(4,2);   %2*pi;
% phiI(8)=x(12)*(xmaxmin(12,1)-xmaxmin(12,2))+xmaxmin(4,2);   %2*pi;

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

% omega=x(3)*10+10;  % [10, 20]  ??
% phi0=x(4)*2*pi;  % [0, 2*pi]
% phiI = zeros(1,8);  % arbitrary phase [0, 2*pi]
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

% after getting M2 and N1, we use 'fzero' function to look for the
% Lagrangian multiplier solution of the nonlinear constraints

B=diag([1,-1,0,0,1,-1,0,0]);
% C is the constraint matrix
C=[1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0; ...
   0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0; ...
   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0; ...
   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0];

lb=-5.0;  % lower boundary for the searching area
ub=5.0;  % upper boundary for the searching area
xzero=[];  % variable contains the roots of function
ftmp=[];
atmp=[];
j=0;
mul=0;  % multiplier

norm1=norm(N1,Inf);
norm2=norm(M2,Inf);

for i=lb:0.02:ub  % here 'i' is not an integer
    if MultiFunc1(i-0.01,M2,N1,norm1,norm2) * MultiFunc1(i+0.01,M2,N1,norm1,norm2) < 0
        % here is one root in the interval
        j=j+1;
        xzero(j)=fzero(@(x)MultiFunc1(x,M2,N1,norm1,norm2),[i-0.01 i+0.01]);
    end 
end

if j>0
    xzero=xzero*norm2;
    for k=1:1:j
        mul=xzero(k);
        atmp=(M-mul*B)\N;
        ftmp(k)=atmp'*N-0.5*atmp'*M*atmp;
    end
end

[fval,I]=max(ftmp);  % here we need maximum value
mul=xzero(I);
a=(M2-mul*B)\N1;


%out1=N1.*(M2\N1)
%out2=M2\N1
%out3=N1

%[LLR,~]=linearConst(N1,M2);
[LLR,~]=nonlinearConst(N1,M2);

% B=diag([1,-1,0,0,1,-1,0,0]);
% % C is the constraint matrix
% C=[1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0; ...
%    0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0; ...
%    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0; ...
%    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0];
% 
% % inverse square matrix, warning message is printed if X is badly scaled or nearly singular.
% a1=M2\N1;  % \ backslash or matrix left division, more efficient
% a2=M2\C';
% a=a1-a2*((C*a2)\(C*a1));
% fbest = -0.5 * ( dot(N1,a1) - (C*a1)'*((C*a2)\(C*a1)) );
% %M2inv=zeros(8,8);
% %M2inv=inv(M2);
% 
% % LLR with 4 constraints, here is accually LLR = - LLR
% if a'*B*a > 0
%    %LLR=-0.5*dot(N1,a1);  % LLR without constraints
%    %LLR = -0.5 * ( N1'*M2inv*N1 - (C*M2inv*N1)'*(C*M2inv*C')*(C*M2inv*N1) );
%    %LLR = -0.5 * ( dot(N1,a1) - (C*a1)'*inv(C*M2inv*C')*(C*a1) );
%    %LLR = -0.5 * ( dot(N1,a1) - (C*a1)'*((C*(M2\C'))\(C*a1)) );  % maybe most efficient and accurate
%    LLR = fbest;
%    %LLR = -0.5 * ( dot(N1,a1) - (C*a1)'*inv(C*a2)*(C*a1) );
%    %LLR = -0.5 * ( dot(N1,a1) - (C*a1)'*(C*a2)\(C*a1) );
% 
% else
%     
%     %disp('Warning: In LogLikelihoodRatio.m: a(1)^2+a(5)^2-a(2)^2-a(6)^2 < 0');
%     LLR=fbest;
% 
% end


% ================

% MM=M2\N1;
% B=diag([1,-1,0,0,1,-1,0,0]);
% 
% a=MM;
% AmpB=0.5*(sqrt(a(1)^2+a(5)^2)+sqrt(a(1)^2+a(5)^2-a(2)^2-a(6)^2));
% 
% if -(a(5)-a(7))/(a(1)-a(3))>0
%     thetaNB=0.5*atan(-(a(5)-a(7))/(a(1)-a(3)));
% elseif -(a(5)-a(7))/(a(1)-a(3))<0
%     thetaNB=0.5*(atan(-(a(5)-a(7))/(a(1)-a(3)))+pi);
% end
% 
% iotacr=-(a(2)+a(4))/2.0/(2.0*AmpB*sin(2.0*thetaNB));
% 
% if MM'*B*MM > 0
%     if iotacr>1.0 || iotacr<-1.0
%         disp('iota may be complex number now !!')
%         LLR=Inf;
%     else
%         LLR=-0.5*dot(N1,MM);  % \ backslash or matrix left division, more efficient
%     end
% else
%    LLR=Inf;
% end

% END of function