% nonlinear constraint multiplier finder function
% Yan Wang, May 27, 2013

function fval=MultiFunc1(mu,M,N)

B=diag([1,-1,0,0,1,-1,0,0]);
% D=diag([0,1,0,0,-1,0,0,0]);

a=(M-mu*B)\N;

fval=a'*B*a;

% END of function