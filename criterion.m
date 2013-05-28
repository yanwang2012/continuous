% integrate the log likelihood ratio, marginalization
% Yan Wang, April 25, 2013.

%clear all;
%clc;

tic;  % Measure performance using stopwatch timer

disp('criterion.m, starting the calculation ... ');

%global Np alphaP deltaP kp N timingResiduals sd yr stdTrueCoord


% set up the parameters and function handles
Nsample=2*10^5;  % number of sampling for each chain
Nchain=2;  % number of Markov chain
distDims=12;
stepsize=0.01;  % stepsize for the proposal distribution
sdprop=0.1;  % standard deviation/Sigma for the proposed Gaussian distribution
start=sdprop*rand(Nchain,distDims); 

% log target distribution
logpdf=@(x) log( mvnpdf(x,zeros(1,12),diag(sdprop*ones(1,12))) );
%logpdf = @(x) LogLikelihoodRatio(x,inParams);

% log proposal distribution, need to be carefully designed
%logproppdf=@(x,y) log(mvnpdf(x-y,[0,0,0],[2,0,0;0,2,0;0,0,2]));
logproppdf=@(x,y) log( mvnpdf(x-y,zeros(1,12),diag(sdprop*ones(1,12))) );

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
%matlabpool open 2

% run the Metropolis-Hastings sampling algorithm
disp('starting the sampling ... ');
% smpl(Nsample,distDims,Nchain), accept counts for the burnin data 
[smpl,accept]=MCMC_MHp(inParamsMH);

%save smpl accept  % save a good realization of the target function for
%future usage

% Shut down the parallel environment
%matlabpool close

% now start to integrate the multi-dimensional function
disp('starting the integration ... ');
burnin=1000;  % burn in period are discarded in the sample
if burnin > Nsample
    disp('Warning: burnin should be much less than Nsample.')
end

% MCMC integration to calculate the volume under the LLR around center
Npoint=6;  % total numer of points to check
intgrt=zeros(Npoint,1);

alphaR=0.2;
deltaR=0.2;
omegaR=0.2;
phi0R=0.2;
phiIR(1:8,1)=0.2;

% create input parameter structure for MCMC_MH()
range = struct('distDims',distDims,'alphaR',alphaR,'deltaR',deltaR,'omegaR',omegaR,'phi0R',phi0R,'phiIR',phiIR);
% initialize the structure, they are the integration range of each
% parameter
range.distDims=distDims;
range.alphaR=alphaR;
range.deltaR=deltaR;
range.omegaR=omegaR;
range.phi0R=phi0R;
range.phiIR(1:8,1)=phiIR(1:8,1);
 
for i=1:1:Npoint
    center=manymins(i).X;
    %center=bestLocationVec(i);
    for l=1:1:Nchain
        
        for j=(burnin+1):1:Nsample
            
            smpltmp=smpl(j,:,l);
            smpltmp=smpltmp+center;  % in standard unit
            
            if rangefunc(smpltmp,range) <= 1.0
                %intgrt(i) = LogLikelihoodRatio(smpltmp,inParams) / exp( logpdf(smpltmp) );
                intgrt(i) = LogLikelihoodRatio(smpltmp,inParams) - logpdf(smpltmp);
            end
            
        end
        
    end
    
    intgrt(i) = intgrt(i)/(Nchain*(Nsample-burnin));
    disp(['integration for ', num2str(i), 'is: ', num2str(intgrt(i)) ]);
    
end

%intgrt = intgrt/(Nchain*(Nsample-burnin));

%disp(num2str(intgrt'));


toc;
% END OF SCRIPT