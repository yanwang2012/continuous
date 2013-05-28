% solve the optimization problem for CW search of PTA by Particle Swarm
% Optimization (PSO)
% April 24, 2013

% notes: in this version, we modified the LoglikelihoodRatio.m, so that we
% are now solving a constraint optimization problem, rather than a
% un-constraint optimization in the previous code. The constraints include
% a(1)+a(3)=0; a(2)-a(4)=0; a(5)+a(7)=0; a(6)-a(8)=0.

%clear all;
%clc;

tic;  % Measure performance using stopwatch timer

disp('testFitness2.m, starting the calculation ... ');

global Np alphaP deltaP kp N timingResiduals sd yr stdTrueCoord

[alpha,delta,omega,phi0,phiI]=PTAsimulator;  % function generating the simulated timing residuals for each pulsar

% set the range of the parameters
xmaxmin=zeros(2,2);  % x_max, x_min for each parameter x, Npara by 2
xmaxmin(1,1)=1.9+0.2; xmaxmin(1,2)=1.9-0.2;
xmaxmin(2,1)=0.6+0.2; xmaxmin(2,2)=0.6-0.2;
% -------
xmaxmin(3,1)=20.0; xmaxmin(3,2)=10.0;
xmaxmin(4,1)=pi; xmaxmin(4,2)=0.0;
xmaxmin(5,1)=pi; xmaxmin(5,2)=0.0;
xmaxmin(6,1)=pi; xmaxmin(6,2)=0.0;
xmaxmin(7,1)=pi; xmaxmin(7,2)=0.0;
xmaxmin(8,1)=pi; xmaxmin(8,2)=0.0;
xmaxmin(9,1)=pi; xmaxmin(9,2)=0.0;
xmaxmin(10,1)=pi; xmaxmin(10,2)=0.0;
xmaxmin(11,1)=pi; xmaxmin(11,2)=0.0;
xmaxmin(12,1)=pi; xmaxmin(12,2)=0.0;

%Standardized true parameter values
stdTrueCoord = zeros(1,12);  % 4+Np 
stdTrueCoord(1)= (alpha-xmaxmin(1,2))/(xmaxmin(1,1)-xmaxmin(1,2));  % [0, 2*pi]
stdTrueCoord(2)= (delta-xmaxmin(2,2))/(xmaxmin(2,1)-xmaxmin(2,2));  % [-pi/2, pi/2]
stdTrueCoord(3)= (omega-10)/10;  % [10, 20]
stdTrueCoord(4)= phi0/pi;  % [0, 2*pi]
stdTrueCoord(5:12)= phiI'/pi;  % [0, 2*pi]


inParams = struct('Np',Np,'N',N,'s',timingResiduals,'sd',sd,...
                  'alphaP',alphaP,'deltaP',deltaP,'kp',kp,'yr',yr,...
                  'xmaxmin',xmaxmin);


% nDims = 12;  % number of dimensions for the parameter space
% xVec = rand(1,nDims);  % randomly choosing initial value for parameters            

% LLR_PSO accepts only the standard coordinates [0,1] as the input variable
fitVal = LLR_PSO(stdTrueCoord,inParams);  % - LogLikelihoodRatio, minimization
disp(['Fitness for perfect subtraction ', num2str(fitVal)]);  % fitVal for true parameters


% -------------------------------
% PARTICLE SWARM OPTIMIZATION
fHandle = @(x) LLR_PSO(x,inParams);

% PSO configuration parameter structure
P = psoparamstruct(1,'default');
P.tuningVars.numPart=20;
P.convergeVars.stepsInCube=30;
P.modScheme.schemeName = '';
P.tuningVars.moveType = 'mv4rand';
defaultModParams = optimset;
defaultModParams.TolX = 0.05;
defaultModParams.TolFun = 0.1;
P.modScheme.schemeParams = defaultModParams;
%%
% PSO output parameter structure
outP = struct('outFileNames',{{'log',''}},'graphics','off','status','off');
%%
% Call PSO
nDim = 12;
nRuns = 10;
bestLocationVec = zeros(nRuns,nDim);
bestFitValVec = zeros(nRuns,1);
for lprun = 1:nRuns
    psoResults=pso(fHandle,nDim,P,outP);  % auto random initial value ?
    disp(['Fitness for run ',num2str(lprun),' is ',num2str(psoResults.bestSNR)]);
    bestLocationVec(lprun,:)=psoResults.bestLocation;
    bestFitValVec(lprun) = psoResults.bestSNR;
end
[bestFitVal,bestFitIndx] = min(bestFitValVec);
disp(['Best Run ',num2str(bestFitIndx)]);
bestLocation = bestLocationVec(bestFitIndx,:);

[~,realC]=LLR_PSO(bestLocation,inParams);
disp('Best Coordinates (before fminsearch)');
disp(realC');
disp(['best fitness found (before fminsearch)', num2str(bestFitVal)]);
% Do fminsearch (Nelder-Mead)
bestLocation = fminsearch(fHandle,bestLocation);
[bestFitVal,realC]=LLR_PSO(bestLocation,inParams);
disp(['best fitness found (after fminsearch)', num2str(bestFitVal)]);
disp('Best Coordinates (after fminsearch)');
disp(realC');

% fitVal for best fit parameters realC

% calculate the maximized extrinsic parameters from the best estimated intrinsic parameters
N1B=zeros(8,1);  % 1-D vector
M2B=zeros(8,8);  % 2-D matrix
AbaseB=zeros(Np,8,N);  % Np*8*N cube, base functions for each pulsar

alphaB=realC(1);
deltaB=realC(2);
omegaB=realC(3);
phi0B=realC(4);

kB=zeros(1,3);  % unit vector pointing from SSB to source
kB(1)=cos(deltaB)*cos(alphaB);
kB(2)=cos(deltaB)*sin(alphaB);
kB(3)=sin(deltaB);

for i=1:1:Np
    AbaseB(i,:,:)=AbaseFunc(alphaB,deltaB,omegaB,phi0B,...
        realC(4+i),alphaP(i),deltaP(i),yr);
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

B=diag([1,-1,0,0,1,-1,0,0]);
% C is the constraint matrix
C=[1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0; ...
   0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0; ...
   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0; ...
   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0];

% inverse square matrix, warning message is printed if X is badly scaled or nearly singular.
a1=M2\N1;  % \ backslash or matrix left division, more efficient
a2=M2\C';
a=a1-a2*((C*a2)\(C*a1));
fbest = -0.5 * ( dot(N1,a1) - (C*a1)'*((C*a2)\(C*a1)) );

if a(1)^2+a(5)^2-a(2)^2-a(6)^2 < 0
   disp('Warning: In testFitness2.m: a(1)^2+a(5)^2-a(2)^2-a(6)^2 < 0');
end
 
AmpB=0.5*(sqrt(a(1)^2+a(5)^2)+sqrt(a(1)^2+a(5)^2-a(2)^2-a(6)^2));

if a(2)^2-4*AmpB^2*a(5)^2/(a(1)^2+a(5)^2) > 0
   disp('Warning: in GOtest1.m: a(2)^2-4*Amp^2*a(5)^2/(a(1)^2+a(5)^2) > 0');
end

if -(a(5)-a(7))/(a(1)-a(3))>0
    thetaNB=0.5*atan(-(a(5)-a(7))/(a(1)-a(3)));
elseif -(a(5)-a(7))/(a(1)-a(3))<0
    thetaNB=0.5*(atan(-(a(5)-a(7))/(a(1)-a(3)))+pi);
end

%iotaB=acos(-a(2)/(2.0*AmpB*sin(2.0*thetaNB)));
iotaB=acos(-(a(2)+a(4))/2.0/(2.0*AmpB*sin(2.0*thetaNB)));

% if -a(5)/a(1)>0
%     thetaNB=0.5*atan(-a(5)/a(1));
% elseif -a(5)/a(1)<0
%     thetaNB=0.5*(atan(-a(5)/a(1))+pi);
% end
% 
% iotaB=acos(-a(2)/(2.0*AmpB*sin(2.0*thetaNB)));

bestfitResiduals=zeros(Np,N);
%is realC a strct or array?? it's 1 by 12 array
for i=1:1:Np
    thetaB=acos(kB*kp(i,:)');
    bestfitResiduals(i,:)=FullResiduals(alphaB,deltaB,omegaB,phi0B,...
    realC(4+i),alphaP(i),deltaP(i),AmpB,iotaB,thetaNB,thetaB,yr);
end

% for i=1:1:Np
%     for mu=1:1:8
%         tmp(1,:)=AbaseB(i,mu,:);
%         bestfitResiduals(i,:)=bestfitResiduals(i,:)+a(mu)*tmp;
%     end
% end

% if a(1)^2+a(5)^2-a(2)^2-a(6)^2>0
%     AmpB=0.5*(sqrt(a(1)^2+a(5)^2)+sqrt(a(1)^2+a(5)^2-a(2)^2-a(6)^2));
% else
%     disp('Warning: a(1)^2+a(5)^2-a(2)^2-a(6)^2 < 0');
%     AmpB=sqrt(a(1)^2+a(5)^2)/(1+cos(iotaB)^2);
% end

% -------------------------------
% plot the timing residuals for each pulsar
figure
for i=1:1:Np
    subplot(4,2,i)
    plot(yr,timingResiduals(i,:),'.-');
    grid on;
    hold on
    plot(yr,bestfitResiduals(i,:),'r.-');
    xlabel('years');
    ylabel('timing residuals (sec)');
    title(pname(i));
    
end

toc;

% END OF SCRIPT