% log likelihood ratio function for full parameter set
% April 18, 2013.  Yan Wang

function [LLR,varargout]=LogLikelihoodRatioF(x,~)

global Np alphaP deltaP kp N timingResiduals sd yr %stdTrueCoord

alpha=x(1)*2*pi;  % [0, 2*pi]
delta=x(2)*pi-pi/2.0;  % [-pi/2, pi/2]
omega=x(3)*10+10;  % [10, 20]  ??
phi0=x(4)*2*pi;  % [0, 2*pi]

phiI = zeros(1,8);  % arbitrary phase [0, 2*pi]
phiI(1)=x(5)*2*pi;
phiI(2)=x(6)*2*pi;
phiI(3)=x(7)*2*pi;
phiI(4)=x(8)*2*pi;
phiI(5)=x(9)*2*pi;
phiI(6)=x(10)*2*pi;
phiI(7)=x(11)*2*pi;
phiI(8)=x(12)*2*pi;

%Amp=10^(x(13));
iota=x(13)*pi;
thetaN=x(14)*pi;

% sky location in Cartesian coordinate
k=zeros(1,3);  % unit vector pointing from SSB to source
k(1)=cos(delta)*cos(alpha);
k(2)=cos(delta)*sin(alpha);
k(3)=sin(delta);

%LLR=0;
tr=zeros(Np,N);

for i=1:1:Np

    theta=acos(k*kp(i,:)');
    %sprintf('%d pulsar theta=%g',i,theta)
    %phiI(i)=mod(phi0-omega*distP(i)*(1-cos(theta)), 2*pi);  % modulus after division
    
    tr(i,:)=FullResiduals(alpha,delta,omega,phi0,phiI(i),alphaP(i),deltaP(i),...
                         1.0,iota,thetaN,theta,yr);  % Amp=1.0

end

tmp1=0;
tmp2=0;
for i=1:1:Np
    tmp1=tmp1+InnProduct(timingResiduals(i,:),tr(i,:),sd);
    tmp2=tmp2+InnProduct(tr(i,:),tr(i,:),sd);
end

Amp=tmp1/tmp2;  % maximized likelihood estimation of Amp

if Amp > 0
    LLR = - 0.5*tmp1^2/tmp2;  % minus sign for the miminization
else
    LLR=Inf;
end


if nargout > 1
    varargout{1}=[alpha,delta,omega,phi0,phiI,iota,thetaN,Amp];
end

%LLR = LLR + InnProduct(timingResiduals(i,:),tr(i,:),sd) ... 
%           -0.5*InnProduct(tr(i,:),tr(i,:),sd);
    