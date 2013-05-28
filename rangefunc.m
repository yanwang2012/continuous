% function to calculate the distance from center to the sampled point
% May 8, 2013

function dist=rangefunc(x,range)

dim=range.distDims;
ra=zeros(dim,1);
ra(1,1)=range.alphaR;
ra(2,1)=range.deltaR;
ra(3,1)=range.omegaR;
ra(4,1)=range.phi0R;
ra(5:dim,1)=range.phiIR;

for i=1:1:dim
    dist=x(1,i)^2/ra(i,1)^2;
end
% END of function