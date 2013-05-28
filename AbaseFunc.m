% Function to calculate A, i.e. the base function
% Yan Wang, March 4, 2013.

%function r=FullResiduals(alpha,delta,omega,phi0,phiI,alphaP,deltaP,Amp,iota,thetaN,theta,t)
function A=AbaseFunc(alpha,delta,omega,phi0,phiI,alphaP,deltaP,t)

N=length(t);
A=zeros(8,N);  % base functions, 8 and N swiched from FullResiduals.m

alphatilde=alpha-alphaP;

Pp=-cos(deltaP)^2*(1-2*cos(alphatilde)^2+cos(alphatilde)^2*cos(delta)^2)...
   +sin(deltaP)^2*cos(delta)^2-0.5*sin(2*deltaP)*cos(alphatilde)*sin(2*delta);

Pc=2*cos(deltaP)*sin(alphatilde)*(cos(deltaP)*cos(alphatilde)*sin(delta)...
   -sin(deltaP)*cos(delta));

p1=[cos(deltaP)*cos(alphaP) cos(deltaP)*sin(alphaP) sin(deltaP)];
p2=[cos(delta)*cos(alpha) cos(delta)*sin(alpha) sin(delta)];
theta=acos(dot(p1,p2));

Fp=Pp/(1-cos(theta));
Fc=Pc/(1-cos(theta));


for i=1:1:N
    
    A(1,i)=Fp*(cos(2*phi0)-cos(2*phiI))*cos(omega*t(i));
    A(2,i)=Fp*(sin(2*phi0)-sin(2*phiI))*cos(omega*t(i));
    A(3,i)=Fp*(sin(2*phi0)-sin(2*phiI))*sin(omega*t(i));
    A(4,i)=Fp*(cos(2*phi0)-cos(2*phiI))*sin(omega*t(i));
    A(5,i)=Fc*(cos(2*phi0)-cos(2*phiI))*cos(omega*t(i));
    A(6,i)=Fc*(sin(2*phi0)-sin(2*phiI))*cos(omega*t(i));
    A(7,i)=Fc*(sin(2*phi0)-sin(2*phiI))*sin(omega*t(i));
    A(8,i)=Fc*(cos(2*phi0)-cos(2*phiI))*sin(omega*t(i));
    
end

% END of function