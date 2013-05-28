% integrate the log likelihood ratio, marginalization
% Yan Wang, April 25, 2013.

clear all;
clc;

tic;  % Measure performance using stopwatch timer

disp('margLLR.m, starting the calculation ... ');

global Np alphaP deltaP kp N timingResiduals sd yr stdTrueCoord

% function generating the simulated timing residuals for each pulsar
[alpha,delta,omega,phi0,phiI]=PTAsimulator;

% set the range of the parameters
xmaxmin=zeros(2,2);  % x_max, x_min for each parameter x, Npara by 2
xmaxmin(1,1)=2*pi; xmaxmin(1,2)=0.0;
xmaxmin(2,1)=pi/2; xmaxmin(2,2)=-pi/2;
% xmaxmin(3,1)=16+3.0; xmaxmin(3,2)=16-5.0;
% xmaxmin(4,1)=1.6+0.2; xmaxmin(4,2)=1.6-0.2;
% xmaxmin(5,1)=3.3+0.2; xmaxmin(5,2)=3.3-0.2;
% xmaxmin(6,1)=4.7+0.2; xmaxmin(6,2)=4.7-0.2;
% xmaxmin(7,1)=0.33+0.2; xmaxmin(7,2)=0.33-0.2;
% xmaxmin(8,1)=4.5+0.2; xmaxmin(8,2)=4.5-0.2;
% xmaxmin(9,1)=0.55+0.2; xmaxmin(9,2)=0.55-0.2;
% xmaxmin(10,1)=1.5+0.2; xmaxmin(10,2)=1.54-0.2;
% xmaxmin(11,1)=0.57+0.2; xmaxmin(11,2)=0.57-0.2;
% xmaxmin(12,1)=0.34+0.2; xmaxmin(12,2)=0.34-0.2;

%Standardized true parameter values
stdTrueCoord = zeros(1,12);  % 4+Np 
stdTrueCoord(1)= (alpha-xmaxmin(1,2))/(xmaxmin(1,1)-xmaxmin(1,2));  % [0, 2*pi]
stdTrueCoord(2)= (delta-xmaxmin(2,2))/(xmaxmin(2,1)-xmaxmin(2,2));  % [-pi/2, pi/2]
stdTrueCoord(3)= (omega-10)/10;  % [10, 20]
stdTrueCoord(4)= phi0/(2*pi);  % [0, 2*pi]
stdTrueCoord(5:12)= phiI'/(2*pi);  % [0, 2*pi]

% create input parameter structure for LogLikelihoodRatio()
inParams = struct('Np',Np,'N',N,'s',timingResiduals,'sd',sd,...
                  'alphaP',alphaP,'deltaP',deltaP,'kp',kp,'yr',yr,'xmaxmin',xmaxmin);

fitVal = LLR_PSO(stdTrueCoord,inParams);  % - LogLikelihoodRatio, minimization
disp(['Fitness for perfect subtraction ', num2str(fitVal)]);  % fitVal for true parameters

% set up the parameters and function handles

Nsample=2*10^5;  % number of sampling for each chain
Nchain=2;  % number of Markov chain
distDims=12;
start=rand(Nchain,distDims);  % starting point, Nchain by distDims
stepsize=0.1;  % stepsize for the proposal distribution

% log target distribution
logpdf=@(x) log( mvnpdf(x,zeros(1,12)+0.5,diag(0.5*ones(1,12))) );
%logpdf = @(x) LogLikelihoodRatio(x,inParams);

% log proposal distribution, need to be carefully designed
%logproppdf=@(x,y) log(mvnpdf(x-y,[0,0,0],[2,0,0;0,2,0;0,0,2]));
logproppdf=@(x,y) log( mvnpdf(x-y,zeros(1,12),diag(1.0*ones(1,12))) );

% propose the random steps
proprnd=@(x,distDims,stepsize) x+rand(1,distDims)*2*stepsize-stepsize;

% create input parameter structure for MCMC_MH()
inParamsMH = struct('start',start,'Nsample',Nsample,'Nchain',Nchain, ...
      'stepsize',stepsize,'logpdf',logpdf,'logproppdf',logproppdf,'proprnd',proprnd);
% initialize the structure
inParamsMH.start=start;
inParamsMH.Nsample=Nsample;
inParamsMH.Nchain=Nchain;
inParamsMH.stepsize=stepsize;
inParamsMH.logpdf=logpdf;
inParamsMH.logproppdf=logproppdf;
inParamsMH.logprnd=proprnd;

% Set up parallel processing, change for to parfor in MCMC_MH()
matlabpool open 2

% run the Metropolis-Hastings sampling algorithm
disp('starting the sampling ... ');
% smpl(Nsample,distDims,Nchain), accept counts for the burnin data 
[smpl,accept]=MCMC_MHp(inParamsMH);

% Shut down the parallel environment
matlabpool close


% ------- display the samples --------  not feasible for distDims > 3
% nc=1;  % chain number
% figure;
% %plot3(smpl(:,1,1),smpl(:,2,1),smpl(:,3,1));
% %plot3(smpl(:,1,nc),smpl(:,2,nc),smpl(:,3,nc),'.r');
% plot(smpl(:,1,nc),'.r');
% grid on


% now start to integrate the multi-dimensional function
disp('starting the integration ... ');
burnin=1000;  % burn in period are discarded in the sample
if burnin > Nsample
    disp('Warning: burnin should be much less than Nsample.')
end

% draw a 2-D grid for the marginalized LLR of alpha and delta, the skymap
% the sky is devided into a 2-D mesh, with value refers the center
% coordinates of each grid
Nalpha=100;  % number of bins
Ndelta=50;  % number of bins, always an even number
skyalpha=zeros(Nalpha,1);
skydelta=zeros(Ndelta,1);
% up to here, alpha and delta have the standard coordinates [0,1]
for i=1:1:Nalpha
    %skyalpha(i)=(i-0.5)*2*pi/Nalpha;
    skyalpha(i)=(i-0.5)/Nalpha;
end
for j=1:1:Ndelta
    %skydelta(j)=-pi/2.0+(j-0.5)*pi/Ndelta;
    skydelta(j)=(j-0.5)/Ndelta;
end

% MCMC integration of the sampled target distribution
intgrt=zeros(Ndelta, Nalpha);  % delta by alpha
ii=0;
for l=1:1:Nchain
    for i=(burnin+1):1:Nsample
        
        if (smpl(i,1,l)<1.0 && smpl(i,1,l)>0.0) && (smpl(i,2,l)<1.0 && smpl(i,2,l)>0.0)
            j=ceil(smpl(i,1,l)/(1.0/Nalpha));  % index for alpha
            k=ceil(smpl(i,2,l)/(1.0/Ndelta));  %+round(Ndelta/2);  % index for delta
            smpltmp=smpl(i,:,l);  % in standard unit
            intgrt(k,j)= LogLikelihoodRatio(smpltmp,inParams) / exp( logpdf(smpltmp) );
            %intgrt=intgrt+(5*mvnpdf(smpl(i,:,l),[1],[2]))/mvnpdf(smpl(i,:,l),[0],[1]);
            ii=ii+1;
        end
        
    end
end

% now transform the standard coordinates to the physical coordinates
skyalpha=(skyalpha-xmaxmin(1,2))/(xmaxmin(1,1)-xmaxmin(1,2));
skydelta=(skydelta-xmaxmin(2,2))/(xmaxmin(2,1)-xmaxmin(2,2));

%intgrt = intgrt/(Nchain*(Nsample-burnin));

% % illustrate the skymap
% figure
% %surf(alphaEs,deltaEs,MLE);
% surf(skydelta,skyalpha,intgrt);  % marginalized log likelihood ratio
% xlabel('declination, \delta');
% ylabel('right ascension, \alpha');
% title('Skymap of GW source, marginalized log likelihood ratio')

%save;  % stores all variables from the current workspace



toc;
% END OF SCRIPT