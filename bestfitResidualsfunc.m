% function calculate the bestfit residuals given the intrinsic
% parameters and others ...
% May 14, 2013

function [bestfitResiduals,varargout]=bestfitResidualsfunc(x,inParams)

%global Np alphaP deltaP kp N timingResiduals sd yr
global kp

xmaxmin = inParams.xmaxmin;
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

% alphaB=realC(1);
% deltaB=realC(2);
% omegaB=realC(3);
% phi0B=realC(4);

k=zeros(1,3);  % unit vector pointing from SSB to source
k(1)=cos(delta)*cos(alpha);
k(2)=cos(delta)*sin(alpha);
k(3)=sin(delta);

for i=1:1:Np
    Abase(i,:,:)=AbaseFunc(alpha,delta,omega,phi0, ...
        phiI(i),alphaP(i),deltaP(i),yr);
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

%[~,a]=linearConst(N1,M2);
[~,a]=nonlinearConst(N1,M2);

% %B=diag([1,-1,0,0,1,-1,0,0]);
% % C is the constraint matrix
% C=[1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0; ...
%    0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0; ...
%    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0; ...
%    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0];
% %a=zeros(8,1);
% a1=M2B\N1B;
% a2=M2B\C';
% a=a1-a2*((C*a2)\(C*a1));

if a(1)^2+a(5)^2-a(2)^2-a(6)^2 < 0
   disp('Warning: in GOtest1.m: a(1)^2+a(5)^2-a(2)^2-a(6)^2 < 0');
end

Amp=0.5*(sqrt(a(1)^2+a(5)^2)+sqrt(a(1)^2+a(5)^2-a(2)^2-a(6)^2));

if a(2)^2-4*Amp^2*a(5)^2/(a(1)^2+a(5)^2) > 0
   disp('Warning: in GOtest1.m: a(2)^2-4*Amp^2*a(5)^2/(a(1)^2+a(5)^2) > 0');
end
 

if -(a(5)-a(7))/(a(1)-a(3))>0
    thetaN=0.5*atan(-(a(5)-a(7))/(a(1)-a(3)));
elseif -(a(5)-a(7))/(a(1)-a(3))<0
    thetaN=0.5*(atan(-(a(5)-a(7))/(a(1)-a(3)))+pi);
end

%iotaB=acos(-a(2)/(2.0*AmpB*sin(2.0*thetaNB)));
iota=acos(-(a(2)+a(4))/2.0/(2.0*Amp*sin(2.0*thetaN)));

bestfitResiduals=zeros(Np,N);
% is realC a strct or array??
for i=1:1:Np
    theta=acos(k*kp(i,:)');
    bestfitResiduals(i,:)=FullResiduals(alpha,delta,omega,phi0,...
    phiI(i),alphaP(i),deltaP(i),Amp,iota,thetaN,theta,yr);
end


% END of function