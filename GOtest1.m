% simulate timing residuals for pta and search for parameters
% April 17, 2013. Yan Wang


% notes:
% 1. parallel processing: matlab pool, 'UseParallel','always'
% Global Search (not support the UseParallel option), Multi Start

tic;  % Measure performance using stopwatch timer.

clear all;
clc;

disp('GOtest1.m, starting the calculation ... ');

global Np alphaP deltaP kp N timingResiduals sd yr stdTrueCoord

% function generating the simulated timing residuals for each pulsar
[alpha,delta,omega,phi0,phiI,pname]=PTAsimulator;

% set the range of the parameters
xmaxmin=zeros(12,2);  % x_max, x_min for each parameter x, Npara by 2
xmaxmin(1,1)=2*pi; xmaxmin(1,2)=0.0;
xmaxmin(2,1)=pi/2; xmaxmin(2,2)=-pi/2;
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
stdTrueCoord(3)= (omega-xmaxmin(3,2))/(xmaxmin(3,1)-xmaxmin(3,2));  % [10, 20]
%stdTrueCoord(4)= phi0/pi;  % [0, 2*pi]
%stdTrueCoord(5:12)= phiI'/pi;  % [0, 2*pi]
stdTrueCoord(4)= mod(phi0,pi)/pi;  % [0, 2*pi]
stdTrueCoord(5:12)= mod(phiI,pi)'/pi;  % [0, 2*pi]

inParams = struct('Np',Np,'N',N,'s',timingResiduals,'sd',sd,...
                  'alphaP',alphaP,'deltaP',deltaP,'kp',kp,'yr',yr,'xmaxmin',xmaxmin);


fitVal = LLR_PSO(stdTrueCoord,inParams);  % - LogLikelihoodRatio, minimization
disp(['Fitness for perfect subtraction ', num2str(fitVal)]);  % fitVal for true parameters

% ---------------------------------------
% GLOBAL OPTIMIZATION

%fHandle = @(x) LLR_PSO(x',inParams);
%[fHandle, varargout] = @(x) LogLikelihoodRatio(x,inParams);
%fHandle = @(x) LogLikelihoodRatio(x,inParams);
fHandle = @(x) LLR_PSO(x,inParams);  % non-linear constraints
%fHandle = @(x) LogLikelihoodRatio(x);

%fHandle = @(x) x(1)^2+x(2)^2+x(3)^2+x(4)^2+x(5)^2+x(6)^2 ... 
%             + x(7)^2+x(8)^2+x(9)^2+x(10)^2+x(11)^2+x(12)^2;

% Create a solver object
ms = MultiStart('Display','iter','TolFun',1e-3,'TolX',1e-3,'MaxTime',6000,...
    'UseParallel','always');
%gs = GlobalSearch('Display','iter','TolFun',1e-3,'TolX',1e-3,'MaxTime',1000);

% Set Start Points for MultiStart (MultiStart only)
%xstart=rand(12,1);
xstart=rand(1,12);
%xstart=stdTrueCoord+0.1*rand(1,12)-0.05;


% Create a problem structure
opts = optimset('Algorithm','interior-point'); % Create an options structure 
%pro=@(x)(2.1*x(1)^2+1.4*x(2)^2);
%pro=@(x)((sin(x(1)-0.1)+sin(x(2)-0.1))/(x(1)^2+x(2)^2));
% fitness=@(x)LLR_PSO(x,inParams);
problem = createOptimProblem('fmincon', ... 
    'x0', xstart, ...
    'objective',fHandle,...
    'lb',[ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. ], ...
    'ub',[ 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. ], ...
    'options',opts);

% Set up parallel processing
matlabpool open 4

% Run the solver
instance=500;  % number of runs
[xmin,fmin,flag,outpt,manymins] = run(ms,problem,instance);
%[xmin,fmin,flag,outpt,manymins] = run(gs,problem);

%sprintf('xmin=%g, fmin=%g, flag=%g, outpt=%g, manymins=%g',xmin,fmin,flag,outpt,manymins)
%xmin,fmin,flag,outpt,manymins

% Shut down the parallel environment
matlabpool close

% for multistart, we could get a set of local minimums from the structure
% 'manymins', so that it's possible that we make the 'bestfit' residuals by
% different minimums we choose
i=1;  % the index for the minimum you want, default is 1
%xminfirst=manymins(i).X;
% calculate the timing residuals for a choice of the bestLocationVec
intPstd=manymins(i).X;  % bestLocation;
%intPstd=bestFitValVec(bestFitIndx,:);  % input intrinsic parameters
bestfitResiduals=bestfitResidualsfunc(intPstd,inParams);

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


% [fbest,a,realC]=LLR_PSO(xminfirst,inParams);
% disp(['Best Coordinates for ', num2str(i), 'th minimum in manymins']);
% disp(realC');
% disp(['best fitness found ', num2str(fbest)]);
% 
% 
% % % calculate the maximized extrinsic parameters from the estimated intrinsic parameters
% % N1B=zeros(8,1);  % 1-D vector
% % M2B=zeros(8,8);  % 2-D matrix
% % AbaseB=zeros(Np,8,N);  % Np*8*N cube, base functions for each pulsar
% % 
% alphaB=realC(1);
% deltaB=realC(2);
% omegaB=realC(3);
% phi0B=realC(4);
% 
% kB=zeros(1,3);  % unit vector pointing from SSB to source
% kB(1)=cos(deltaB)*cos(alphaB);
% kB(2)=cos(deltaB)*sin(alphaB);
% kB(3)=sin(deltaB);
% % 
% % for i=1:1:Np
% %     AbaseB(i,:,:)=AbaseFunc(alphaB,deltaB,omegaB,phi0B,...
% %         realC(4+i),alphaP(i),deltaP(i),yr);
% % end
% % 
% % tmp=zeros(1,N);
% % for mu=1:1:8
% % 
% %     % calculate 8*1 vector N
% %     for i=1:1:Np
% %         tmp(1,:)=AbaseB(i,mu,:);  % convert (1,1,N) matrix to (1,N)
% %         %size(tmp)
% %         N1B(mu)=N1B(mu)+InnProduct(timingResiduals(i,:),tmp,sd);
% %         %N1(mu)=N1(mu)+InnProduct(tmp1,tmp2,sd);
% %     end
% %     
% %     % calculate 8*8 matrix M
% %     for nu=1:1:8
% %         
% %         for i=1:1:Np
% %             M2B(mu,nu)=M2B(mu,nu)+InnProduct(AbaseB(i,mu,:),AbaseB(i,nu,:),sd);
% %         end
% %         
% %     end
% % 
% % end
% % 
% % %B=diag([1,-1,0,0,1,-1,0,0]);
% % % C is the constraint matrix
% % C=[1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0; ...
% %    0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0; ...
% %    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0; ...
% %    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0];
% % %a=zeros(8,1);
% % a1=M2B\N1B;
% % a2=M2B\C';
% % a=a1-a2*((C*a2)\(C*a1));
% 
% if a(1)^2+a(5)^2-a(2)^2-a(6)^2 < 0
%    disp('Warning: in GOtest1.m: a(1)^2+a(5)^2-a(2)^2-a(6)^2 < 0');
% end
% 
% AmpB=0.5*(sqrt(a(1)^2+a(5)^2)+sqrt(a(1)^2+a(5)^2-a(2)^2-a(6)^2));
% 
% if a(2)^2-4*AmpB^2*a(5)^2/(a(1)^2+a(5)^2) > 0
%    disp('Warning: in GOtest1.m: a(2)^2-4*Amp^2*a(5)^2/(a(1)^2+a(5)^2) > 0');
% end
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
% 
% bestfitResiduals=zeros(Np,N);
% % is realC a strct or array??
% for i=1:1:Np
%     thetaB=acos(kB*kp(i,:)');
%     bestfitResiduals(i,:)=FullResiduals(alphaB,deltaB,omegaB,phi0B,...
%     realC(4+i),alphaP(i),deltaP(i),AmpB,iotaB,thetaNB,thetaB,yr);
% end
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
%     xlabel('years');
%     ylabel('timing residuals (sec)');
%     title(pname(i));
%     
% end

% Display the results and relevant statistics
% figure,
% fcnvals=[manymins.Fval];
% hist([manymins.Fval]);
% grid on
% ylim([0 5])


toc;
% END OF SCRIPT