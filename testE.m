% test for extrinsic parameter estimations

global Np alphaP deltaP kp N timingResiduals sd yr % stdTrueCoord

PTAsimulator;

% calculate the maximized extrinsic parameters from the estimated intrinsic parameters
N1B=zeros(8,1);  % 1-D vector
M2B=zeros(8,8);  % 2-D matrix
AbaseB=zeros(Np,8,N);  % Np*8*N cube, base functions for each pulsar

alphaB=pi/4;
deltaB=pi/4;
omegaB=2*pi/0.3925;
phi0B=0.8;

kB=zeros(1,3);  % unit vector pointing from SSB to source
kB(1)=cos(deltaB)*cos(alphaB);
kB(2)=cos(deltaB)*sin(alphaB);
kB(3)=sin(deltaB);

phiI=zeros(8,1);
phiI(1,1)=6.1269;
phiI(2,1)=4.0791;
phiI(3,1)=3.1626;
phiI(4,1)=0.2953;
phiI(5,1)=0.0071;
phiI(6,1)=4.6230;
phiI(7,1)=5.6629;
phiI(8,1)=4.6521;


for i=1:1:Np
    AbaseB(i,:,:)=AbaseFunc(alphaB,deltaB,omegaB,phi0B,...
        phiI(i),alphaP(i),deltaP(i),yr);
end

tmp=zeros(1,N);
for mu=1:1:8

    % calculate 8*1 vector N
    for i=1:1:Np
        tmp(1,:)=AbaseB(i,mu,:);  % convert (1,1,N) matrix to (1,N)
        %size(tmp)
        N1B(mu)=N1B(mu)+InnProduct(timingResiduals(i,:),tmp,sd);
        %N1(mu)=N1(mu)+InnProduct(tmp1,tmp2,sd);
    end
    
    % calculate 8*8 matrix M
    for nu=1:1:8
        
        for i=1:1:Np
            M2B(mu,nu)=M2B(mu,nu)+InnProduct(AbaseB(i,mu,:),AbaseB(i,nu,:),sd);
        end
        
    end

end

%a=zeros(8,1);
a=M2B\N1B;

AmpB=0.5*(sqrt(a(1)^2+a(5)^2)+sqrt(a(1)^2+a(5)^2-a(2)^2-a(6)^2));

if -a(5)/a(1)>0
    thetaNB=0.5*atan(-a(5)/a(1));
elseif -a(5)/a(1)<0
    thetaNB=0.5*(atan(-a(5)/a(1))+pi);
end

iotaB=acos(-a(2)/(2.0*AmpB*sin(2.0*thetaNB)));

bestfitResiduals=zeros(Np,N);
% % is realC a strct or array??
% for i=1:1:Np
%     thetaB=acos(kB*kp(i,:)');
%     bestfitResiduals(i,:)=FullResiduals(alphaB,deltaB,omegaB,phi0B,...
%     realC(4+i),alphaP(i),deltaP(i),AmpB,iotaB,thetaNB,thetaB,yr);
% end
% 
% 
% % -------------------------------
% % plot the timing residuals for each pulsar
% figure
% for i=1:1:Np
%     subplot(4,2,i)
%     plot(yr,timingResiduals(i,:),'.-');
%     grid on;
%     hold on
%     plot(yr,bestfitResiduals(i,:),'r.-');
%     xlabel('Modified Juliant Day');
%     ylabel('timing residuals (sec)');
%     title('simulated residuals (blue) and bestfit residuals (red)')
%     
% end

