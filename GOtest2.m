% simulate timing residuals for pta and search for parameters
% April 17, 2013. Yan Wang
% use the log likelihood function with full set of parameters

% notes:
% 1. parallel processing: matlab pool, 'UseParallel','always'
% Global Search (not support the UseParallel option), Multi Start

tic;  % Measure performance using stopwatch timer.

clear all;
clc;

disp('GOtest2.m, starting the calculation ... ');

global Np alphaP deltaP kp N timingResiduals sd yr stdTrueCoord

PTAsimulator;  % function generating the simulated timing residuals for each pulsar

inParams = struct('Np',Np,'N',N,'s',timingResiduals,'sd',sd,...
                  'alphaP',alphaP,'deltaP',deltaP,'kp',kp,'yr',yr);

%fitVal = LLR_PSO(stdTrueCoord,inParams);  % - LogLikelihoodRatio, minimization
xinput=zeros(1,14);  % 12+3-1(Amp)
for i=1:1:12
    xinput(i)=stdTrueCoord(i);
end
%xinput(13)=log(1.3*10^(-7));
xinput(13)=(pi/4+0.6)/pi;  % iota
xinput(14)=(pi/4+0.3)/pi;  % thetaN

fitVal = LogLikelihoodRatioF(xinput,inParams);
disp(['Fitness for perfect subtraction ', num2str(fitVal)]);  % fitVal for true parameters
             
% ---------------------------------------
% GLOBAL OPTIMIZATION

fHandle = @(x) LLR_PSO(x,inParams);
%[fHandle, varargout] = @(x) LogLikelihoodRatio(x,inParams);
%fHandle = @(x) LogLikelihoodRatio(x,inParams);
%fHandle = @(x) LogLikelihoodRatioF(x);

%fHandle = @(x) x(1)^2+x(2)^2+x(3)^2+x(4)^2+x(5)^2+x(6)^2 ... 
%             + x(7)^2+x(8)^2+x(9)^2+x(10)^2+x(11)^2+x(12)^2;

% Create a solver object
ms = MultiStart('Display','iter','TolFun',1e-3,'TolX',1e-3,'MaxTime',2000,...
    'UseParallel','always');
%gs = GlobalSearch('Display','iter','TolFun',1e-3,'TolX',1e-3,'MaxTime',1000);

% Set Start Points for MultiStart (MultiStart only)
%xstart=rand(12,1);
%xstart=stdTrueCoord;  %+0.1*rand(1,12)-0.05;
xstart=xinput+0.05*rand(1,14)-0.025;

% Create a problem structure
opts = optimset('Algorithm','interior-point'); % Create an options structure 
%pro=@(x)(2.1*x(1)^2+1.4*x(2)^2);
%pro=@(x)((sin(x(1)-0.1)+sin(x(2)-0.1))/(x(1)^2+x(2)^2));
% fitness=@(x)LLR_PSO(x,inParams);
problem = createOptimProblem('fmincon', ... 
    'x0', xstart, ...
    'objective',fHandle,...
    'lb',[ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. ], ...
    'ub',[ 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. ], ...
    'options',opts);

% Set up parallel processing
%matlabpool open 4

% Run the solver
[xmin,fmin,flag,outpt,manymins] = run(ms,problem,100);
%[xmin,fmin,flag,outpt,manymins] = run(gs,problem);

%sprintf('xmin=%g, fmin=%g, flag=%g, outpt=%g, manymins=%g',xmin,fmin,flag,outpt,manymins)
%xmin,fmin,flag,outpt,manymins

% Shut down the parallel environment
%matlabpool close

%[fbest,realC]=LLR_PSO(xmin,inParams);
[fbest,realC]=LogLikelihoodRatioF(xmin,inParams);
disp('Best Coordinates');
disp(realC');
disp(['best fitness found ', num2str(fbest)]);


% calculate the maximized extrinsic parameters from the estimated intrinsic parameters
% N1B=zeros(8,1);  % 1-D vector
% M2B=zeros(8,8);  % 2-D matrix
% AbaseB=zeros(Np,8,N);  % Np*8*N cube, base functions for each pulsar

alphaB=realC(1);
deltaB=realC(2);
omegaB=realC(3);
phi0B=realC(4);

kB=zeros(1,3);  % unit vector pointing from SSB to source
kB(1)=cos(deltaB)*cos(alphaB);
kB(2)=cos(deltaB)*sin(alphaB);
kB(3)=sin(deltaB);

% for i=1:1:Np
%     AbaseB(i,:,:)=AbaseFunc(alphaB,deltaB,omegaB,phi0B,...
%         realC(4+i),alphaP(i),deltaP(i),yr);
% end
% 
% tmp=zeros(1,N);
% for mu=1:1:8
% 
%     % calculate 8*1 vector N
%     for i=1:1:Np
%         tmp(1,:)=AbaseB(i,mu,:);  % convert (1,1,N) matrix to (1,N)
%         %size(tmp)
%         N1B(mu)=N1B(mu)+InnProduct(timingResiduals(i,:),tmp,sd);
%         %N1(mu)=N1(mu)+InnProduct(tmp1,tmp2,sd);
%     end
%     
%     % calculate 8*8 matrix M
%     for nu=1:1:8
%         
%         for i=1:1:Np
%             M2B(mu,nu)=M2B(mu,nu)+InnProduct(AbaseB(i,mu,:),AbaseB(i,nu,:),sd);
%         end
%         
%     end
% 
% end
% 
% %a=zeros(8,1);
% a=M2B\N1B;
% 
% if a(1)^2+a(5)^2-a(2)^2-a(6)^2 < 0
%    disp('Warning: a(1)^2+a(5)^2-a(2)^2-a(6)^2 < 0');
% end
%  
% AmpB=0.5*(sqrt(a(1)^2+a(5)^2)+sqrt(a(1)^2+a(5)^2-a(2)^2-a(6)^2));
% 
%     
% if -(a(5)-a(7))/(a(1)-a(3))>0
%     thetaNB=0.5*atan(-(a(5)-a(7))/(a(1)-a(3)));
% elseif -(a(5)-a(7))/(a(1)-a(3))<0
%     thetaNB=0.5*(atan(-(a(5)-a(7))/(a(1)-a(3)))+pi);
% end
% 
% %iotaB=acos(-a(2)/(2.0*AmpB*sin(2.0*thetaNB)));
% iotaB=acos(-(a(2)+a(4))/2.0/(2.0*AmpB*sin(2.0*thetaNB)));

bestfitResiduals=zeros(Np,N);
% is realC a strct or array??
for i=1:1:Np
    thetaB=acos(kB*kp(i,:)');
    bestfitResiduals(i,:)=FullResiduals(alphaB,deltaB,omegaB,phi0B,...
    realC(4+i),alphaP(i),deltaP(i),realC(15),realC(13),realC(14),thetaB,yr);
end

% -------------------------------
% plot the timing residuals for each pulsar
figure
for i=1:1:Np
    subplot(4,2,i)
    plot(yr,timingResiduals(i,:),'.-');
    grid on;
    hold on
    plot(yr,bestfitResiduals(i,:),'r.-');
    xlabel('Modified Juliant Day');
    ylabel('timing residuals (sec)');
    title('simulated residuals (blue) and bestfit residuals (red)')
    
end

% Display the results and relevant statistics
% figure,
% fcnvals=[manymins.Fval];
% hist([manymins.Fval]);
% grid on
% ylim([0 5])


toc;
% END OF SCRIPT