% integration/marginalization of a multi-dimensional parameter functions
% within certain ranges of the variables
% Yan Wang, April 20, 2013.

clear all;
clc;

tic;  % Measure performance using stopwatch timer.

disp('integf2.m, starting the calculation ... ');

% set up the parameters and function handles
start=[1.1; -3.4];  % starting point, Nchain by Nsample
Nsample=2*10^5;  % number of sampling
Nchain=size(start,1);  % number of Markov chain
stepsize=0.5;  % stepsize for the proposal distribution

%distDims=size(start,2);

% log target distribution
%logpdf=@(x) log( mvnpdf(x,[0,0,0],[1,0,0;0,1,0;0,0,1])...
%        +mvnpdf(x,[5,5,5],[2,0,0;0,2,0;0,0,2]) );
logpdf=@(x) log( mvnpdf(x,[0],[1]) );

% log proposal distribution, need to be carefully designed
%logproppdf=@(x,y) log(mvnpdf(x-y,[0,0,0],[2,0,0;0,2,0;0,0,2]));
logproppdf=@(x,y) log(mvnpdf(x-y,[0],[2]));

% propose the random steps
proprnd=@(x,distDims,stepsize) x+rand(1,distDims)*2*stepsize-stepsize;

% create input parameter structure for MCMC_MH()
inParams = struct('start',start,'Nsample',Nsample,'Nchain',Nchain, ...
      'stepsize',stepsize,'logpdf',logpdf,'logproppdf',logproppdf,'proprnd',proprnd);
% initialize the structure
inParams.start=start;
inParams.Nsample=Nsample;
inParams.Nchain=Nchain;
inParams.stepsize=stepsize;
inParams.logpdf=logpdf;
inParams.logproppdf=logproppdf;
inParams.logprnd=proprnd;

% Set up parallel processing, change for to parfor in MCMC_MH()
matlabpool open 2

% run the Metropolis-Hastings sampling algorithm
% smpl(Nsample,distDims,Nchain), accept counts for the burnin data 
[smpl,accept]=MCMC_MHp(inParams);

% Shut down the parallel environment
matlabpool close


% ------- display the samples --------
nc=1;  % chain number
figure;
%plot3(smpl(:,1,1),smpl(:,2,1),smpl(:,3,1));
%plot3(smpl(:,1,nc),smpl(:,2,nc),smpl(:,3,nc),'.r');
plot(smpl(:,1,nc),'.r');
grid on

% now start to integrate the multi-dimensional function
burnin=1000;
if burnin > Nsample
    disp('Warning: burnin should be much less than Nsample.')
end

% MCMC integration of the sampled target distribution
intgrt=0;
ii=0;
for l=1:1:Nchain
    for i=(burnin+1):1:Nsample
        
        intgrt=intgrt+(5*mvnpdf(smpl(i,:,l),[1],[2]))/mvnpdf(smpl(i,:,l),[0],[1]);
        ii=ii+1;
        
    end
end

intgrt = intgrt/(Nchain*(Nsample-burnin));

    
toc;
% END OF SCRIPT

% notes:
% MHSAMPLE implements Independent Metropolis-Hasting sampling.