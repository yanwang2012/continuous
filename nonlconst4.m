% non-linear constraint functions
% May 1, 2013

function [c, ceq] = nonlconst4(a)

% Nonlinear inequality constraints, c(a)<=0
c = [];
%c = [-a(1)^2-a(5)^2+a(2)^2+a(6)^2;
%    a(2)^2-(sqrt(a(1)^2+a(5)^2)+sqrt(a(1)^2+a(5)^2-a(2)^2-a(6)^2))^2*sin(2*thetaN)^2];
% c = [ -a(1)^2-a(5)^2+a(2)^2+a(6)^2;
%         a(2)^2-(sqrt(a(1)^2+a(5)^2)+sqrt(a(1)^2+a(5)^2-a(2)^2-a(6)^2))^2 ...
%         *a(5)^2/(a(1)^2+a(5)^2) ];

%c = -a(1)^2-a(5)^2+a(2)^2+a(6)^2;

%c = -a(1)^2-a(3)^2+a(2)^2+a(4)^2;

% c = [1.5 + x(1)*x(2) - x(1) - x(2);     
%      -x(1)*x(2) - 10];
 
% Nonlinear equality constraints
%ceq = [ -a(1)^2-a(5)^2+a(2)^2+a(6)^2;
%        a(2)^2-(sqrt(a(1)^2+a(5)^2)+sqrt(a(1)^2+a(5)^2-a(2)^2-a(6)^2))^2 ...
%        *a(5)^2/(a(1)^2+a(5)^2) ];
    
% ceq = [ -a(1)^2-a(5)^2+a(2)^2+a(6)^2;
%         a(2)^2-a(5)^2 ];

%ceq = a(2)^2-a(5)^2;
ceq = [ a(1)*a(2)+a(3)*a(4);
        -a(1)^2-a(3)^2+a(2)^2+a(4)^2];
%ceq = [];
% END OF FUNCTION