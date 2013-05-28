% non-linear constraint functions
% April 29, 2013

function [c, ceq] = nonlconst1(a)

% Amp=0.5*(sqrt(a(1)^2+a(5)^2)+sqrt(a(1)^2+a(5)^2-a(2)^2-a(6)^2));

% if -(a(5)-a(7))/(a(1)-a(3))>0
%     thetaN=0.5*atan(-(a(5)-a(7))/(a(1)-a(3)));
% elseif -(a(5)-a(7))/(a(1)-a(3))<0
%     thetaN=0.5*(atan(-(a(5)-a(7))/(a(1)-a(3)))+pi);
% end

% Nonlinear inequality constraints, c(a)<=0
%c = [];
%c = -a(1)^2-a(5)^2+a(2)^2+a(6)^2;

c = -a(1)^2-a(3)^2+a(2)^2+a(4)^2;
%c = [-a(1)^2-a(5)^2+a(2)^2+a(6)^2;
%     a(2)^2-(sqrt(a(1)^2+a(5)^2)+sqrt(a(1)^2+a(5)^2-a(2)^2-a(6)^2))^2*sin(2*thetaN)^2];
%c = -a(1)^2-a(5)^2+a(2)^2+a(6)^2;

% c = [1.5 + x(1)*x(2) - x(1) - x(2);     
%      -x(1)*x(2) - 10];
 
% Nonlinear equality constraints
%ceq = -a(1)^2-a(5)^2+a(2)^2+a(6)^2;
%ceq = -a(1)^2-a(3)^2+a(2)^2+a(4)^2;
ceq=[];

% END OF FUNCTION