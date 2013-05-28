% Function for Markov chain Monte Carlo (Metropolis-Hastings Algorithm)
% This is a basic Metropolis-Hastings sampling algorithm, which samples
% a target multi-dimensional function/distribution with proposal 
% distribution, and return the samples for different chains for future 
% integration or optimization...
% -----------------------------------------------
% Input:
% 'start': the starting point of sampling
% 'inParams' structure
% '.Nsample': the number of samplings
% '.Nchain': the number of chains for each run
% '.stepsize': stepsize for proposal function
% -----------------------------------------------
% Output:
% 'smpl': Nsample by Nchain by distDims matrix
% 'accept': the acception rate for each Markov chain, it measures the 
%           efficiency of the algorithm
%
% -----------------------------------------------
% Yan Wang, March 5, 2013.
% modified: April 19, 2013 - add 'inParams' structure, 'Nchain', parallel

function [smpl,accept]=MCMC_MH(inParams)

% pass variables in 'inParams' structure to local variables
start=inParams.start;
Nsample=inParams.Nsample;
Nchain=inParams.Nchain;
stepsize=inParams.stepsize;
logpdf=inParams.logpdf;
logproppdf=inParams.logproppdf;
proprnd=inParams.proprnd;

outclass=superiorfloat(start);  % single or double         
distDims=size(start,2);  % dimension of the distribution
smpl=zeros([Nsample,distDims,Nchain],outclass);  % output chains

x0=start;  % x0 is the place holder for the current value, initialized to start
y=zeros(Nchain,distDims);

accept=zeros(Nchain,1,outclass);  % accept rate for each chain, 
%sprintf('Nsample=%d, Nchain=%d, distDims=%d',Nsample,Nchain,distDims);
%sprintf('x0=%f',x0);


% Metropolis-Hastings sampling, use logarithm formalism to avoid overflow
u=log(rand(Nsample,Nchain));
q1=zeros(Nchain,1);
q2=zeros(Nchain,1);
rho=zeros(Nchain,1);
r=zeros(Nchain,1);


parfor l=1:1:Nchain  % parallel chains
%for l=1:1:Nchain
    tmp=zeros(Nsample,distDims);
    for i=1:1:Nsample
        y(l,:)=proprnd(x0(l,:),distDims,stepsize);  % propose a new position with stepsize
        q1(l)=logproppdf(x0(l,:),y(l,:));
        q2(l)=logproppdf(y(l,:),x0(l,:));
        
        % Metropolis ratio
        rho(l)=(q1(l)+logpdf(y(l,:)))-(q2(l)+logpdf(x0(l,:)));  % generic formula, nonsymmetric.
        %rho=logpdf(y)-logpdf(x0);  % for symmetric proposal pdf, q1=q2
        
        % accept or reject the proposal
        r(l)=min(rho(l),0);
        %ui=u(i,l);
        
        if r(l)>=u(i,l)
            % accept, move to new position y
            tmp(i,:)=y(l,:);
            x0(l,:)=y(l,:);
            accept(l)=accept(l)+1;
        else
            % reject, stay at current position x0
            tmp(i,:)=y(l,:);
            %smpl(i,l,:)=x0(l,:);
        end
        
    end
    
    smpl(i,:,l)=tmp(i,:);
    
end

% Accept rate can be used to optimize the choice of scale parameters in
% random walk MH sampler. See for example Roberts, Gelman and Gilks (1997).
accept=accept/Nsample;

% Move the replicates dimension to the end to make samples easier to
% manipulate.
smpl = permute(smpl,[1 3 2]);


%{

%-------------------------------------------------
% target pdf function to be sampled
function y=logpdf(x)
  %y=log( normpdf(x,0,1)+normpdf(x,10,3) );
  y=log( mvnpdf(x,[0,0,0],[1,0,0;0,1,0;0,0,1])...
        +mvnpdf(x,[5,5,5],[2,0,0;0,2,0;0,0,2]) );


%-------------------------------------------------
% log proposal distribution, need to be carefully designed.
% In general, it is possible to use suboptimal inference and learning algorithms 
% to generate data-driven proposal distributions. (Andrieu et al, 2003)
function output=logproppdf(x,y)
  %output=log(normpdf(x-y));
  output=log(mvnpdf(x-y,[0,0,0],[2,0,0;0,2,0;0,0,2]));
  

%-------------------------------------------------
% propose the random steps
function y=proprnd(x,distDims,stepsize)
  % stepsize=0.5;
  y=x+rand(1,distDims)*2*stepsize-stepsize;

% [smpl,acc]=MCMC_MH(1.1,10^5,1);
% histfit(smpl,100)
% histo2D(D,[0 5],10,[0 10],5, 'Xlabel','Ylabel','2d Histogram') ;
% plot(smpl(:,1,1),smpl(:,1,2))

%}

% END of function