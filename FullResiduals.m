% Function to calculate the timing residuals for a full set of parameters 
% (initial phase, arbitrary phase ...)
% Yan Wang, March 4, 2013.

function r=FullResiduals(alpha,delta,omega,phi0,phiI,alphaP,deltaP,Amp,iota,thetaN,theta,t)

N=length(t);

r=zeros(N,1);
c=zeros(8,1);  % coefficient function
A=zeros(N,8);  % base functions

%sprintf('len N=%d',N)
%theta=acos(k.*kp);
%theta=acos(dot(k,kp));
%Delta=distP*(1-cos(theta));

alphatilde=alpha-alphaP;

Pp=-cos(deltaP)^2*(1-2*cos(alphatilde)^2+cos(alphatilde)^2*cos(delta)^2)...
   +sin(deltaP)^2*cos(delta)^2-0.5*sin(2*deltaP)*cos(alphatilde)*sin(2*delta);

Pc=2*cos(deltaP)*sin(alphatilde)*(cos(deltaP)*cos(alphatilde)*sin(delta)...
   -sin(deltaP)*cos(delta));

Fp=Pp/(1-cos(theta));
Fc=Pc/(1-cos(theta));

% here Amp=\zeta*omega^(-1/3) that defined in my notes
c(1,1)=-Amp*(1+cos(iota)^2)*cos(2*thetaN);
c(2,1)=-Amp*2*cos(iota)*sin(2*thetaN);
c(3,1)= Amp*(1+cos(iota)^2)*cos(2*thetaN);
c(4,1)=-Amp*2*cos(iota)*sin(2*thetaN);
c(5,1)= Amp*(1+cos(iota)^2)*sin(2*thetaN);
c(6,1)=-Amp*2*cos(iota)*cos(2*thetaN);
c(7,1)=-Amp*(1+cos(iota)^2)*sin(2*thetaN);    
c(8,1)=-Amp*2*cos(iota)*cos(2*thetaN);

%phiI=phi0-omega*distP*(1-cos(theta));

for i=1:1:N
    
    A(i,1)=Fp*(cos(2*phi0)-cos(2*phiI))*cos(omega*t(i));
    A(i,2)=Fp*(sin(2*phi0)-sin(2*phiI))*cos(omega*t(i));
    A(i,3)=Fp*(sin(2*phi0)-sin(2*phiI))*sin(omega*t(i));
    A(i,4)=Fp*(cos(2*phi0)-cos(2*phiI))*sin(omega*t(i));
    A(i,5)=Fc*(cos(2*phi0)-cos(2*phiI))*cos(omega*t(i));
    A(i,6)=Fc*(sin(2*phi0)-sin(2*phiI))*cos(omega*t(i));
    A(i,7)=Fc*(sin(2*phi0)-sin(2*phiI))*sin(omega*t(i));
    A(i,8)=Fc*(cos(2*phi0)-cos(2*phiI))*sin(omega*t(i));
    
    for j=1:1:8 
        r(i)=r(i)+c(j)*A(i,j);
    end
    
end

% END of function